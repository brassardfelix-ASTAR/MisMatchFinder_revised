# MismatchFinder – Split Source (Inline Restructure)

This folder contains the code split from `CHATGPT.MMFdebug.rs` into a Rust module tree.

## Layout
- `src/main.rs`
- `src/lib.rs`
- `src/bamreader.rs` (top-level module; declares submodules below)
- `src/bamreader/mismatch.rs`
- `src/bamreader/fragments.rs`
- `src/bamreader/cigar.rs`
- `src/bamreader/filter.rs` (declares `germline` and `region`)
- `src/bamreader/filter/germline.rs`
- `src/bamreader/filter/region.rs`
- `src/output.rs` (declares `tsvwriter` and `vcfwriter`)
- `src/output/vcfwriter.rs`

⚠️ **Note:** The original aggregate did not include `output/tsvwriter.rs`. The module is referenced by `src/output.rs` but missing in source. You can:
- Add your existing `tsvwriter.rs` file under `src/output/`, or
- Temporarily stub it:

```rust
// src/output/tsvwriter.rs
pub fn write_mismatches<T>(_x: T) -> Result<(), Box<dyn std::error::Error>> {
    unimplemented!("tsvwriter module was not included in the provided aggregate file");
}
```
