// Copyright 2014 Christopher Schröder, Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.


pub mod record;
pub mod header;
pub mod pileup;
pub mod buffer;

use std::ffi;
use std::ptr;
use std::slice;
use std::path::Path;
use url::Url;
use libc;

use htslib;

pub use bam::record::Record;
pub use bam::header::Header;
pub use bam::buffer::RecordBuffer;


/// A trait for a BAM reader with a read method.
pub trait Read: Sized {
    /// Read next BAM record into given record.
    /// Use this method in combination with a single allocated record to avoid the reallocations
    /// occurring with the iterator.
    ///
    /// # Arguments
    ///
    /// * `record` - the record to be filled
    fn read(&mut self, record: &mut record::Record) -> Result<(), ReadError>;

    /// Iterator over the records of the seeked region.
    /// Note that, while being convenient, this is less efficient than pre-allocating a
    /// `Record` and reading into it with the `read` method, since every iteration involves
    /// the allocation of a new `Record`.
    fn records(&mut self) -> Records<Self>;

    /// Iterator over pileups.
    fn pileup(&mut self) -> pileup::Pileups<Self>;

    /// Return the BGZF struct
    fn bgzf(&self) -> *mut htslib::Struct_BGZF;

    /// Return the header.
    fn header(&self) -> &HeaderView;

    /// Seek to the given virtual offset in the file
    fn seek(&mut self, offset: i64) -> Result<(), SeekError> {
        let ret = unsafe { htslib::bgzf_seek(self.bgzf(), offset, libc::SEEK_SET) };
        if ret == 0 {
             Ok(())
        } else {
            Err(SeekError::Some)
        }
    }

    /// Report the current virtual offset
    fn tell(&self) -> i64  {
        unsafe { htslib::bgzf_tell_func(self.bgzf()) }
    }
}


/// A BAM reader.
pub struct Reader {
    bgzf: *mut htslib::Struct_BGZF,
    header: HeaderView,
}


unsafe impl Send for Reader {}


impl Reader {
    /// Create a new Reader from path.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, ReaderPathError> {
        match path.as_ref().to_str() {
            Some(p) if path.as_ref().exists() => {
                Ok(try!(Self::new(p.as_bytes())))
            },
            _ => {
                Err(ReaderPathError::InvalidPath)
            }
        }
    }

    /// Create a new Reader from STDIN.
    pub fn from_stdin() -> Result<Self, BGZFError> {
        Self::new(b"-")
    }

    /// Create a new Reader from URL.
    pub fn from_url(url: &Url) -> Result<Self, BGZFError> {
        Self::new(url.as_str().as_bytes())
    }

    /// Create a new Reader.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open. Use "-" for stdin.
    fn new(path: &[u8]) -> Result<Self, BGZFError> {
        let bgzf = try!(bgzf_open(&ffi::CString::new(path).unwrap(), b"r"));
        let header = unsafe { htslib::bam_hdr_read(bgzf) };
        Ok(Reader { bgzf: bgzf, header: HeaderView::new(header) })
    }

    extern fn pileup_read(data: *mut ::libc::c_void, record: *mut htslib::bam1_t) -> ::libc::c_int {
        let _self = unsafe { &*(data as *mut Self) };
        unsafe { htslib::bam_read1(_self.bgzf, record) }
    }
}


impl Read for Reader {
    fn read(&mut self, record: &mut record::Record) -> Result<(), ReadError> {
        match unsafe { htslib::bam_read1(self.bgzf, record.inner) } {
            -1 => Err(ReadError::NoMoreRecord),
            -2 => Err(ReadError::Truncated),
            -4 => Err(ReadError::Invalid),
            _  => Ok(())
        }
    }

    /// Iterator over the records of the fetched region.
    /// Note that, while being convenient, this is less efficient than pre-allocating a
    /// `Record` and reading into it with the `read` method, since every iteration involves
    /// the allocation of a new `Record`.
    fn records(&mut self) -> Records<Self> {
        Records { reader: self }
    }

    fn pileup(&mut self) -> pileup::Pileups<Self> {
        let _self = self as *const Self;
        let itr = unsafe {
            htslib::bam_plp_init(
                Some(Reader::pileup_read),
                _self as *mut ::libc::c_void
            )
        };
        pileup::Pileups::new(self, itr)
    }

    fn bgzf(&self) -> *mut htslib::Struct_BGZF {
        self.bgzf
    }

    fn header(&self) -> &HeaderView {
        &self.header
    }
}


impl Drop for Reader {
    fn drop(&mut self) {
        unsafe {
            htslib::bgzf_close(self.bgzf);
        }
    }
}


pub struct IndexedReader {
    bgzf: *mut htslib::Struct_BGZF,
    header: HeaderView,
    idx: *mut htslib::hts_idx_t,
    itr: Option<*mut htslib:: hts_itr_t>,
}


unsafe impl Send for IndexedReader {}


impl IndexedReader {

    /// Create a new Reader from path.
    ///
    /// # Arguments
    ///
    /// * `path` - the path to open.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, IndexedReaderPathError> {
        match path.as_ref().to_str() {
            Some(p) if path.as_ref().exists() => {
                Ok(try!(Self::new(&ffi::CString::new(p).unwrap())))
            },
            _ => {
                Err(IndexedReaderPathError::InvalidPath)
            }
        }
    }

    pub fn from_url(url: &Url) -> Result<Self, IndexedReaderError> {
        Self::new(&ffi::CString::new(url.as_str()).unwrap())
    }

    /// Create a new Reader.
    ///
    /// # Arguments
    ///
    /// * `path` - the path. Use "-" for stdin.
    fn new(path: &ffi::CStr) -> Result<Self, IndexedReaderError> {
        let bgzf = try!(bgzf_open(path, b"r"));
        let header = unsafe { htslib::bam_hdr_read(bgzf) };
        let idx = unsafe {
            htslib::hts_idx_load(
                path.as_ptr(),
                htslib::HTS_FMT_BAI
            )
        };
        if idx.is_null() {
            Err(IndexedReaderError::InvalidIndex)
        } else {
            Ok(IndexedReader { bgzf: bgzf, header: HeaderView::new(header), idx: idx, itr: None })
        }
    }

    pub fn fetch(&mut self, tid: u32, beg: u32, end: u32) -> Result<(), FetchError> {
        if let Some(itr) = self.itr {
            unsafe { htslib::hts_itr_destroy(itr) }
        }
        let itr = unsafe {
            htslib::sam_itr_queryi(self.idx, tid as i32, beg as i32, end as i32)
        };
        if itr.is_null() {
            self.itr = None;
            Err(FetchError::Some)
        }
        else {
            self.itr = Some(itr);
            Ok(())
        }
    }

    extern fn pileup_read(data: *mut ::libc::c_void, record: *mut htslib::bam1_t) -> ::libc::c_int {
        let _self = unsafe { &*(data as *mut Self) };
        match _self.itr {
            Some(itr) => itr_next(_self.bgzf, itr, record),
            None      => 0
        }
    }
}


impl Read for IndexedReader {
    fn read(&mut self, record: &mut record::Record) -> Result<(), ReadError> {
        match self.itr {
            Some(itr) => match itr_next(self.bgzf, itr, record.inner) {
                -1 => Err(ReadError::NoMoreRecord),
                -2 => Err(ReadError::Truncated),
                -4 => Err(ReadError::Invalid),
                _  => Ok(())
            },
            None      => Err(ReadError::NoMoreRecord)
        }
    }

    /// Iterator over the records of the fetched region.
    /// Note that, while being convenient, this is less efficient than pre-allocating a
    /// `Record` and reading into it with the `read` method, since every iteration involves
    /// the allocation of a new `Record`.
    fn records(&mut self) -> Records<Self> {
        Records { reader: self }
    }

    fn pileup(&mut self) -> pileup::Pileups<Self> {
        let _self = self as *const Self;
        let itr = unsafe {
            htslib::bam_plp_init(
                Some(IndexedReader::pileup_read),
                _self as *mut ::libc::c_void
            )
        };
        pileup::Pileups::new(self, itr)
    }

    fn bgzf(&self) -> *mut htslib::Struct_BGZF {
        self.bgzf
    }

    fn header(&self) -> &HeaderView {
        &self.header
    }
}


impl Drop for IndexedReader {
    fn drop(&mut self) {
        unsafe {
            if self.itr.is_some() {
                htslib::hts_itr_destroy(self.itr.unwrap());
            }
            htslib::hts_idx_destroy(self.idx);
            htslib::bgzf_close(self.bgzf);
        }
    }
}


/// A BAM writer.
pub struct Writer {
    f: *mut htslib::Struct_BGZF,
    header: HeaderView,
}


unsafe impl Send for Writer {}


impl Writer {
    /// Create a new BAM file.
    ///
    /// # Arguments
    ///
    /// * `path` - the path.
    /// * `header` - header definition to use
    pub fn from_path<P: AsRef<Path>>(path: P, header: &header::Header) -> Result<Self, WriterPathError> {
        if let Some(p) = path.as_ref().to_str() {
                Ok(try!(Self::new(p.as_bytes(), b"w", header)))
        } else {
            Err(WriterPathError::InvalidPath)
        }
    }

    /// Create a new BAM file at STDOUT.
    ///
    /// # Arguments
    ///
    /// * `header` - header definition to use
    pub fn from_stdout(header: &header::Header) -> Result<Self, BGZFError> {
        Self::new(b"-", b"w", header)
    }

    /// Create a new BAM file.
    ///
    /// # Arguments
    ///
    /// * `path` - the path. Use "-" for stdin.
    /// * `header` - header definition to use
    pub fn new(path: &[u8], mode: &[u8], header: &header::Header) -> Result<Self, BGZFError> {
        let f = try!(bgzf_open(&ffi::CString::new(path).unwrap(), mode));

        // sam_hdr_parse does not populate the text and l_text fields of the header_record.
        // This causes non-SQ headers to be dropped in the output BAM file.
        // To avoid this, we copy the full header to a new C-string that is allocated with malloc,
        // and set this into header_record manually.
        let header_record = unsafe {
            let mut header_string = header.to_bytes();
            if !header_string.is_empty() && header_string[header_string.len() - 1] != b'\n' {
                header_string.push(b'\n');
            }
            let l_text = header_string.len();
            let text = ::libc::malloc(l_text + 1);
            libc::memset(text, 0, l_text + 1);
            libc::memcpy(text, header_string.as_ptr() as *const ::libc::c_void, header_string.len());

            //println!("{}", str::from_utf8(&header_string).unwrap());
            let rec = htslib::sam_hdr_parse(
                (l_text + 1) as i32,
                text as *const i8,
            );

            (*rec).text = text as *mut i8;
            (*rec).l_text = l_text as u32;
            rec
        };

        unsafe { htslib::bam_hdr_write(f, header_record); }

        Ok(Writer { f: f, header: HeaderView::new(header_record) })
    }

    /// Activate multi-threaded BAM write support in htslib. This should permit faster
    /// writing of large BAM files.
    /// # Arguments
    ///
    /// * `n_threads` - number of background writer threads to use
    pub fn set_threads(&mut self, n_threads: usize) -> Result<(), ThreadingError> {
        let r = unsafe { htslib::bgzf_mt(self.f, n_threads as ::libc::c_int, 256) };
        if r != 0 {
            Err(ThreadingError::Some)
        } else {
            Ok(())
        }
    }

    /// Write record to BAM.
    ///
    /// # Arguments
    ///
    /// * `record` - the record to write
    pub fn write(&mut self, record: &record::Record) -> Result<(), WriteError> {
        if unsafe { htslib::bam_write1(self.f, record.inner) } == -1 {
            Err(WriteError::Some)
        }
        else {
            Ok(())
        }
    }


    /// Return the header.
    pub fn header(&self) -> &HeaderView {
        &self.header
    }
}


impl Drop for Writer {
    fn drop(&mut self) {
        unsafe {
            htslib::bgzf_close(self.f);
        }
    }
}


/// Iterator over the records of a BAM.
pub struct Records<'a, R: 'a + Read> {
    reader: &'a mut R
}


impl<'a, R: Read> Iterator for Records<'a, R> {
    type Item = Result<record::Record, ReadError>;

    fn next(&mut self) -> Option<Result<record::Record, ReadError>> {
        let mut record = record::Record::new();
        match self.reader.read(&mut record) {
            Err(ReadError::NoMoreRecord) => None,
            Ok(())   => Some(Ok(record)),
            Err(err) => Some(Err(err))
        }
    }
}


quick_error! {
    #[derive(Debug)]
    pub enum ReadError {
        Truncated {
            description("truncated record")
        }
        Invalid {
            description("invalid record")
        }
        NoMoreRecord {
            description("no more record")
        }
    }
}


impl ReadError {
    /// Returns true if no record has been read because the end of the file was reached.
    pub fn is_eof(&self) -> bool {
        match self {
            &ReadError::NoMoreRecord => true,
            _ => false
        }
    }
}


quick_error! {
    #[derive(Debug)]
    pub enum IndexedReaderError {
        InvalidIndex {
            description("invalid index")
        }
        BGZFError(err: BGZFError) {
            from()
        }
    }
}


quick_error! {
    #[derive(Debug)]
    pub enum WriterPathError {
        InvalidPath {
            description("invalid path")
        }
        BGZFError(err: BGZFError) {
            from()
        }
    }
}


quick_error! {
    #[derive(Debug)]
    pub enum IndexedReaderPathError {
        InvalidPath {
            description("invalid path")
        }
        IndexedReaderError(err: IndexedReaderError) {
            from()
        }
    }
}


quick_error! {
    #[derive(Debug)]
    pub enum BGZFError {
        Some {
            description("error reading BGZF file")
        }
    }
}


quick_error! {
    #[derive(Debug)]
    pub enum ReaderPathError {
        InvalidPath {
            description("invalid path")
        }
        BGZFError(err: BGZFError) {
            from()
        }
    }
}


quick_error! {
    #[derive(Debug)]
    pub enum ThreadingError {
        Some {
            description("error setting threads for multi-threaded writing")
        }
    }
}

quick_error! {
    #[derive(Debug)]
    pub enum WriteError {
        Some {
            description("error writing record")
        }
    }
}


quick_error! {
    #[derive(Debug)]
    pub enum FetchError {
        Some {
            description("error fetching a locus")
        }
    }
}

quick_error! {
    #[derive(Debug)]
    pub enum SeekError {
        Some {
            description("error seeking to voffset")
        }
    }
}


/// Wrapper for opening a BAM file.
fn bgzf_open(path: &ffi::CStr, mode: &[u8]) -> Result<*mut htslib::Struct_BGZF, BGZFError> {
    let ret = unsafe {
        htslib::bgzf_open(
            path.as_ptr(),
            ffi::CString::new(mode).unwrap().as_ptr()
        )
    };
    if ret.is_null() {
        Err(BGZFError::Some)
    } else {
        Ok(ret)
    }
}


/// Wrapper for iterating an indexed BAM file.
fn itr_next(bgzf: *mut htslib::Struct_BGZF, itr: *mut htslib:: hts_itr_t, record: *mut htslib::bam1_t) -> i32 {
    unsafe {
        htslib::hts_itr_next(
            bgzf,
            itr,
            record as *mut ::libc::c_void,
            ptr::null_mut()
        )
    }
}


pub struct HeaderView {
    inner: *mut htslib::bam_hdr_t,
    owned: bool,
}


impl HeaderView {
    /// Create a new HeaderView from the underlying Htslib type, and own it.
    pub fn new(inner: *mut htslib::bam_hdr_t) -> Self {
        HeaderView {
            inner: inner,
            owned: true,
        }
    }

    #[inline]
    pub fn inner(&self) -> htslib::bam_hdr_t {
        unsafe { (*self.inner) }
    }

    pub fn tid(&self, name: &[u8]) -> Option<u32> {
        let tid = unsafe {
            htslib::bam_name2id(
                self.inner,
                ffi::CString::new(name).ok().expect("Expected valid name.").as_ptr()
            )
        };
        if tid < 0 {
            None
        }
        else {
            Some(tid as u32)
        }
    }

    pub fn target_count(&self) -> u32 {
        self.inner().n_targets as u32
    }

    pub fn target_names(&self) -> Vec<&[u8]> {
        let names = unsafe { slice::from_raw_parts(self.inner().target_name, self.target_count() as usize) };
        names.iter().map(|name| unsafe { ffi::CStr::from_ptr(*name).to_bytes() }).collect()
    }

    pub fn target_len(&self, tid: u32) -> Option<u32> {
        let inner = unsafe { *self.inner };
        if (tid as i32) < inner.n_targets {
            let l: &[u32] = unsafe { slice::from_raw_parts(inner.target_len, inner.n_targets as usize) };
            Some(l[tid as usize])
        }
        else {
            None
        }
    }

    /// Retrieve the textual SAM header as bytes
    pub fn as_bytes<'a>(&'a self) -> &'a [u8] {
        unsafe{ ffi::CStr::from_ptr((*self.inner).text).to_bytes() }
    }
}


impl Clone for HeaderView {
    fn clone(&self) -> Self {
        HeaderView {
            inner: unsafe { htslib::bam_hdr_dup(self.inner) },
            owned: true
        }
    }
}


impl Drop for HeaderView {
    fn drop(&mut self) {
        if self.owned {
            unsafe { htslib::bam_hdr_destroy(self.inner); }
        }
    }
}


#[cfg(test)]
mod tests {
    extern crate tempdir;
    use std::collections::HashMap;
    use super::*;
    use super::record::{Cigar, CigarString, Aux};
    use super::header::HeaderRecord;
    use std::str;
    use std::path::Path;

    fn gold() -> ([&'static [u8]; 6], [u16; 6], [&'static [u8]; 6], [&'static [u8]; 6], [CigarString; 6]) {
        let names = [&b"I"[..], &b"II.14978392"[..], &b"III"[..], &b"IV"[..], &b"V"[..], &b"VI"[..]];
        let flags = [16u16, 16u16, 16u16, 16u16, 16u16, 2048u16];
        let seqs = [
            &b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA"[..],
            &b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA"[..],
            &b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA"[..],
            &b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA"[..],
            &b"CCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA"[..],
            &b"ACTAAGCCTAAGCCTAAGCCTAAGCCAATTATCGATTTCTGAAAAAATTATCGAATTTTCTAGAAATTTTGCAAATTTTTTCATAAAATTATCGATTTTA"[..],
        ];
        let quals = [
            &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..],
            &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..],
            &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..],
            &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..],
            &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..],
            &b"#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..],
        ];
        let cigars = [
            CigarString(vec![Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)]),
            CigarString(vec![Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)]),
            CigarString(vec![Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)]),
            CigarString(vec![Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)]),
            CigarString(vec![Cigar::Match(27), Cigar::Del(1), Cigar::Match(73)]),
            CigarString(vec![Cigar::Match(27), Cigar::Del(100000), Cigar::Match(73)]),
        ];
        (names, flags, seqs, quals, cigars)
    }



    #[test]
    fn test_read() {
        let (names, flags, seqs, quals, cigars) = gold();
        let mut bam = Reader::from_path(&Path::new("test/test.bam")).ok().expect("Error opening file.");
        let del_len = [1, 1, 1, 1, 1, 100000];

        for (i, record) in bam.records().enumerate() {
            let rec = record.ok().expect("Expected valid record");
            println!("{}", str::from_utf8(rec.qname()).ok().unwrap());
            assert_eq!(rec.qname(), names[i]);
            assert_eq!(rec.flags(), flags[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);
            assert_eq!(*rec.cigar(), cigars[i]);
            let cigar = rec.cigar();
            if let Ok(end_pos) = cigar.end_pos() {
                assert_eq!( end_pos, rec.pos() + 100 + del_len[i]);
                assert_eq!(cigar.read_pos( end_pos as u32 - 10, false, false).unwrap().unwrap(), 90);
            } else {
                panic!("bug: failed to fetch cigar.end_pos() in test_read()")
            }
            assert_eq!(cigar.read_pos(rec.pos() as u32 + 20, false, false).unwrap().unwrap(), 20);
            assert_eq!(cigar.read_pos(4000000, false, false).unwrap(), None);
            // fix qual offset
            let qual: Vec<u8> = quals[i].iter().map(|&q| q - 33).collect();
            assert_eq!(rec.qual(), &qual[..]);
        }
    }

    #[test]
    fn test_seek() {
        let mut bam = Reader::from_path(&Path::new("test/test.bam")).ok().expect("Error opening file.");

        let mut names_by_voffset = HashMap::new();

        let mut offset = bam.tell();
        let mut rec = Record::new();
        loop {
            if let Err(e) = bam.read(&mut rec) {
                if e.is_eof() {
                    break;
                } else {
                    panic!("error reading bam");
                }
            }
            let qname = str::from_utf8(rec.qname()).unwrap().to_string();
            println!("{} {}", offset, qname);
            names_by_voffset.insert(offset, qname);
            offset = bam.tell();
        }

        for (offset, qname) in names_by_voffset.iter() {
            println!("{} {}", offset, qname);
            bam.seek(*offset).unwrap();
            bam.read(&mut rec).unwrap();
            let rec_qname = str::from_utf8(rec.qname()).ok().unwrap().to_string();
            assert_eq!(qname, &rec_qname);
        }
    }


    #[test]
    fn test_read_sam_header() {
        let mut bam = Reader::from_path(&"test/test.bam").ok().expect("Error opening file.");

        let true_header = "@SQ\tSN:CHROMOSOME_I\tLN:15072423\n@SQ\tSN:CHROMOSOME_II\tLN:15279345\n@SQ\tSN:CHROMOSOME_III\tLN:13783700\n@SQ\tSN:CHROMOSOME_IV\tLN:17493793\n@SQ\tSN:CHROMOSOME_V\tLN:20924149\n".to_string();
        let header_text = String::from_utf8(bam.header.as_bytes().to_owned()).unwrap();
        assert_eq!(header_text, true_header);
    }


    #[test]
    fn test_read_indexed() {
        let (names, flags, seqs, quals, cigars) = gold();
        let mut bam = IndexedReader::from_path(&"test/test.bam").ok().expect("Expected valid index.");

        let tid = bam.header.tid(b"CHROMOSOME_I").expect("Expected tid.");
        assert!(bam.header.target_len(tid).expect("Expected target len.") == 15072423);

        // fetch to position containing reads
        bam.fetch(tid, 0, 2).ok().expect("Expected successful fetch.");
        assert!(bam.records().count() == 6);

        // compare reads
        bam.fetch(tid, 0, 2).ok().expect("Expected successful fetch.");
        for (i, record) in bam.records().enumerate() {
            let rec = record.ok().expect("Expected valid record");
            println!("{}", str::from_utf8(rec.qname()).ok().unwrap());
            assert_eq!(rec.qname(), names[i]);
            assert_eq!(rec.flags(), flags[i]);
            assert_eq!(rec.seq().as_bytes(), seqs[i]);
            assert_eq!(*rec.cigar(), cigars[i]);
            // fix qual offset
            let qual: Vec<u8> = quals[i].iter().map(|&q| q - 33).collect();
            assert_eq!(rec.qual(), &qual[..]);
            assert_eq!(rec.aux(b"NotAvailableAux"), None);
        }

        // fetch to empty position
        bam.fetch(2, 1, 1).ok().expect("Expected successful fetch.");
        assert!(bam.records().count() == 0);
    }

    #[test]
    fn test_set_record() {

        let (names, _, seqs, quals, cigars) = gold();

        let mut rec = record::Record::new();
        rec.set_reverse();
        rec.set(names[0], &cigars[0], seqs[0], quals[0]);
        // note: this segfaults if you push_aux() before set()
        //       because set() obliterates aux
        rec.push_aux(b"NM", &Aux::Integer(15));

        assert_eq!(rec.qname(), names[0]);
        assert_eq!(*rec.cigar(), cigars[0]);
        assert_eq!(rec.seq().as_bytes(), seqs[0]);
        assert_eq!(rec.qual(), quals[0]);
        assert!(rec.is_reverse());
        assert_eq!(rec.aux(b"NM").unwrap(), Aux::Integer(15));
    }

    #[test]
    fn test_set_qname() {

        let (names, _, seqs, quals, cigars) = gold();

        assert!(names[0] != names[1]);

        let mut rec = record::Record::new();
        rec.set(names[0], &cigars[0], seqs[0], quals[0]);
        rec.push_aux(b"NM", &Aux::Integer(15));

        assert_eq!(rec.qname(), names[0]);
        assert_eq!(*rec.cigar(), cigars[0]);
        assert_eq!(rec.seq().as_bytes(), seqs[0]);
        assert_eq!(rec.qual(), quals[0]);
        assert_eq!(rec.aux(b"NM").unwrap(), Aux::Integer(15));

        // Equal length qname
        assert!(rec.qname()[0] != 'X' as u8);
        rec.set_qname(b"X");
        assert_eq!(rec.qname(), b"X");

        // Longer qname
        let mut longer_name = names[0].to_owned().clone();
        let extension = b"BuffaloBUffaloBUFFaloBUFFAloBUFFALoBUFFALO";
        longer_name.extend(extension.iter());
        rec.set_qname(&longer_name);

        assert_eq!(rec.qname(), longer_name.as_slice());
        assert_eq!(*rec.cigar(), cigars[0]);
        assert_eq!(rec.seq().as_bytes(), seqs[0]);
        assert_eq!(rec.qual(), quals[0]);
        assert_eq!(rec.aux(b"NM").unwrap(), Aux::Integer(15));

        // Shorter qname
        let shorter_name = b"42";
        rec.set_qname(shorter_name);

        assert_eq!(rec.qname(), shorter_name);
        assert_eq!(*rec.cigar(), cigars[0]);
        assert_eq!(rec.seq().as_bytes(), seqs[0]);
        assert_eq!(rec.qual(), quals[0]);
        assert_eq!(rec.aux(b"NM").unwrap(), Aux::Integer(15));

        // Zero-length qname
        rec.set_qname(b"");

        assert_eq!(rec.qname(), b"");
        assert_eq!(*rec.cigar(), cigars[0]);
        assert_eq!(rec.seq().as_bytes(), seqs[0]);
        assert_eq!(rec.qual(), quals[0]);
        assert_eq!(rec.aux(b"NM").unwrap(), Aux::Integer(15));
    }

    #[test]
    fn test_remove_aux() {
        let mut bam = Reader::from_path(&Path::new("test/test.bam")).ok().expect("Error opening file.");

        for record in bam.records() {
            let rec = record.ok().expect("Expected valid record");

            if rec.aux(b"XS").is_some() {
                rec.remove_aux(b"XS");
            }

            if rec.aux(b"YT").is_some() {
	            rec.remove_aux(b"YT");
            }

	    rec.remove_aux(b"ab");

	    assert_eq!(rec.aux(b"XS"), None);
	    assert_eq!(rec.aux(b"YT"), None);
        }
    }



    #[test]
    fn test_write() {
        let (names, _, seqs, quals, cigars) = gold();

        let tmp = tempdir::TempDir::new("rust-htslib").ok().expect("Cannot create temp dir");
        let bampath = tmp.path().join("test.bam");
        println!("{:?}", bampath);
        {
            let mut bam = Writer::from_path(
                &bampath,
                Header::new().push_record(
                    HeaderRecord::new(b"SQ").push_tag(b"SN", &"chr1")
                                            .push_tag(b"LN", &15072423)
                )
            ).ok().expect("Error opening file.");

            for i in 0..names.len() {
                let mut rec = record::Record::new();
                rec.set(names[i], &cigars[i], seqs[i], quals[i]);
                rec.push_aux(b"NM", &Aux::Integer(15));

                bam.write(&mut rec).ok().expect("Failed to write record.");
            }
        }

        {
            let mut bam = Reader::from_path(&bampath).ok().expect("Error opening file.");

            for i in 0..names.len() {
                let mut rec = record::Record::new();
                bam.read(&mut rec).ok().expect("Failed to read record.");

                assert_eq!(rec.qname(), names[i]);
                assert_eq!(*rec.cigar(), cigars[i]);
                assert_eq!(rec.seq().as_bytes(), seqs[i]);
                assert_eq!(rec.qual(), quals[i]);
                assert_eq!(rec.aux(b"NM").unwrap(), Aux::Integer(15));
            }
        }

        tmp.close().ok().expect("Failed to delete temp dir");
    }


    #[test]
    fn test_write_threaded() {
        let (names, _, seqs, quals, cigars) = gold();

        let tmp = tempdir::TempDir::new("rust-htslib").ok().expect("Cannot create temp dir");
        let bampath = tmp.path().join("test.bam");
        println!("{:?}", bampath);
        {
            let mut bam = Writer::from_path(
                &bampath,
                Header::new().push_record(
                    HeaderRecord::new(b"SQ").push_tag(b"SN", &"chr1")
                                            .push_tag(b"LN", &15072423)
                )
            ).ok().expect("Error opening file.");
            bam.set_threads(4).unwrap();

            for i in 0 .. 10000 {
                let mut rec = record::Record::new();
                let idx = i % names.len();
                rec.set(names[idx], &cigars[idx], seqs[idx], quals[idx]);
                rec.push_aux(b"NM", &Aux::Integer(15));
                rec.set_pos(i as i32);

                bam.write(&mut rec).ok().expect("Failed to write record.");
            }
        }

        {
            let mut bam = Reader::from_path(&bampath).ok().expect("Error opening file.");

            for (i, _rec) in bam.records().enumerate() {
                let idx = i % names.len();

                let rec = _rec.expect("Failed to read record.");

                assert_eq!(rec.pos(), i as i32);
                assert_eq!(rec.qname(), names[idx]);
                assert_eq!(*rec.cigar(), cigars[idx]);
                assert_eq!(rec.seq().as_bytes(), seqs[idx]);
                assert_eq!(rec.qual(), quals[idx]);
                assert_eq!(rec.aux(b"NM").unwrap(), Aux::Integer(15));
            }
        }

        tmp.close().ok().expect("Failed to delete temp dir");
    }



    #[test]
    fn test_copy_template() {
        // Verify that BAM headers are transmitted correctly when using an existing BAM as a template for headers.

        let tmp = tempdir::TempDir::new("rust-htslib").ok().expect("Cannot create temp dir");
        let bampath = tmp.path().join("test.bam");
        println!("{:?}", bampath);

        let mut input_bam = Reader::from_path(&"test/test.bam").ok().expect("Error opening file.");

        {
            let mut bam = Writer::from_path(
                &bampath,
                &Header::from_template(&input_bam.header()),
            ).ok().expect("Error opening file.");

            for rec in input_bam.records() {
                bam.write(&rec.unwrap()).ok().expect("Failed to write record.");
            }
        }

        {
            let mut copy_bam = Reader::from_path(&bampath).ok().expect("Error opening file.");

            // Verify that the header came across correctly
            assert_eq!(input_bam.header().as_bytes(), copy_bam.header().as_bytes());
        }

        tmp.close().ok().expect("Failed to delete temp dir");
    }



    #[test]
    fn test_pileup() {
        let (_, _, seqs, quals, _) = gold();

        let mut bam = Reader::from_path(&"test/test.bam").ok().expect("Error opening file.");
        let pileups = bam.pileup();
        for pileup in pileups.take(26) {
            let _pileup = pileup.ok().expect("Expected successful pileup.");
            let pos = _pileup.pos() as usize;
            assert_eq!(_pileup.depth(), 6);
            assert!(_pileup.tid() == 0);
            for (i, a) in _pileup.alignments().enumerate() {
                assert_eq!(a.indel(), pileup::Indel::None);
                let qpos = a.qpos().unwrap();
                assert_eq!(qpos, pos - 1);
                assert_eq!(a.record().seq()[qpos], seqs[i][qpos]);
                assert_eq!(a.record().qual()[qpos], quals[i][qpos] - 33);
            }
        }
    }
}
