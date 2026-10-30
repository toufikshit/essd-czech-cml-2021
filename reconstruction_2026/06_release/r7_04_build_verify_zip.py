"""Stage 7: build ESSD_CZ_CML_2021_dataset.zip from the package folder and verify it.

Layout as the earlier archive: primary/*, supplementary_wvd/*, README_DATA.md, DATA_DICTIONARY.md, SHA256SUMS.txt;
stored without recompression (ZIP64). Verification reopens the archive: member list, sizes, CRC (testzip) and the
SHA-256 of every member read from the archive against SHA256SUMS.txt and against the files on disk.
Usage: python -I r7_04_build_verify_zip.py <package dir> <zip path>
"""
import hashlib, os, sys, zipfile

pkg, zpath = sys.argv[1], sys.argv[2]
members = [f'{sub}/{fn}' for sub in ('primary', 'supplementary_wvd') for fn in sorted(os.listdir(os.path.join(pkg, sub)))]
members += ['README_DATA.md', 'DATA_DICTIONARY.md', 'SHA256SUMS.txt']
if os.path.exists(zpath):
    sys.exit(f'{zpath} exists; not overwritten')
with zipfile.ZipFile(zpath, 'w', compression=zipfile.ZIP_STORED, allowZip64=True, strict_timestamps=False) as z:
    for m in members:
        z.write(os.path.join(pkg, m), m)


def sha_stream(f):
    h = hashlib.sha256()
    for b in iter(lambda: f.read(1 << 24), b''):
        h.update(b)
    return h.hexdigest()


sums = dict(reversed(l.split('  ', 1)) for l in open(os.path.join(pkg, 'SHA256SUMS.txt'), encoding='utf-8').read().split('\n') if l)
ok = True
with zipfile.ZipFile(zpath) as z:
    names = [i.filename for i in z.infolist()]
    ok &= names == members
    print('members in order:', names == members, len(names))
    bad = z.testzip(); ok &= bad is None
    print('CRC test:', 'all members OK' if bad is None else f'FAILED at {bad}')
    for i in z.infolist():
        with z.open(i) as f:
            h = sha_stream(f)
        disk = os.path.getsize(os.path.join(pkg, i.filename))
        ref = sums.get(i.filename)
        good = i.file_size == disk and i.compress_type == zipfile.ZIP_STORED and (ref is None or h == ref)
        ok &= good
        print(f"{'OK ' if good else 'BAD'} {i.file_size:>12} {h} {i.filename}" + ('' if ref else '  (checksum file itself)'))
size = os.path.getsize(zpath)
with open(zpath, 'rb') as f:
    zsha = sha_stream(f)
md5 = hashlib.md5()
with open(zpath, 'rb') as f:
    for b in iter(lambda: f.read(1 << 24), b''):
        md5.update(b)
print(f'ZIP {zpath}\n  size {size:,} bytes\n  sha256 {zsha}\n  md5 {md5.hexdigest()}')
print('ALL PASS' if ok else 'FAILED')
