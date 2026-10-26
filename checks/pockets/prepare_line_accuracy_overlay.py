"""Generate a scoped scalar-Laplacian hook without modifying vendor source."""
import hashlib
import json
from pathlib import Path


def prepare_overlay(vendor,output):
    vendor=Path(vendor);output=Path(output)
    relative=Path('CosmoInterface/evolvers/kernels/scalarsingletkernels.h')
    source=vendor/'include'/relative;text=source.read_text()
    include='#include "CosmoInterface/definitions/potential.h"'
    needle='LatLapl(model.fldS(n))'
    if text.count(needle)!=1 or text.count(include)!=1:
        raise ValueError('Unsupported upstream scalar kernel; no overlay generated')
    changed=text.replace(include,include+'\n#include "line_laplacian.h"').replace(needle,'PocketPlanarSpatialLaplacian(model, n)')
    target=output/relative;target.parent.mkdir(parents=True,exist_ok=False);target.write_text(changed)
    report=dict(vendor_source=str(source.resolve()),vendor_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        overlay=str(target.resolve()),overlay_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        change='Only scalar Laplacian dispatch; potential, PV2 and other sectors remain upstream',
        source_untouched=source.read_text()==text)
    (output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('vendor',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();print(json.dumps(prepare_overlay(a.vendor,a.output),indent=2))
