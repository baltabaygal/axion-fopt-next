"""Generate a narrowly pinned PV timing overlay without modifying upstream."""
import hashlib
from pathlib import Path
import sys
import subprocess

source, output = map(Path, sys.argv[1:])
commit = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
if commit != "0dad16ec2c894a7c6abbea0e98734126b3470b36":
    raise SystemExit("Require the official CosmoLattice v2.0.0 commit")
if subprocess.run(["git", "-C", str(source), "diff", "--quiet", "HEAD", "--",
                   "include", "source", "models", "CMakeLists.txt", "cmake"]).returncode:
    raise SystemExit("CosmoLattice source has tracked modifications; use a clean checkout")
relative = Path("CosmoInterface/evolvers/positionverlet.h")
raw = (source / "include" / relative).read_bytes()
expected = "4ce54a76b83647164b03243e2b38121cbe7fca2963a60546bbb9b97613334ebe"
if hashlib.sha256(raw).hexdigest() != expected:
    raise SystemExit("Unsupported positionverlet.h: require unmodified CosmoLattice v2.0.0")
text = raw.decode()

def insert_before(marker, addition):
    global text
    if text.count(marker) != 1:
        raise SystemExit("Ambiguous upstream hook location")
    text = text.replace(marker, addition + marker)

insert_before("        // Now we compute the kicks (pi_0 --> pi_1):", """        // axion-fopt-next: external background at the force midpoint.
        if constexpr (requires { model.setExternalTime(tcurrent); })
          model.setExternalTime(tcurrent + w * model.dt / 2.0);

""")
insert_before("        if (expansion && !fixedBackground) storeFieldsAveragesInt(model);", """        // axion-fopt-next: lazily refresh the endpoint mass if it is measured.
        if constexpr (requires { model.deferExternalTime(tcurrent); })
          model.deferExternalTime(tcurrent + w * model.dt);
        else if constexpr (requires { model.setExternalTime(tcurrent); })
          model.setExternalTime(tcurrent + w * model.dt);

""")
insert_before("      if (fixedBackground) model.aDotI = aBackground.dot(tMinust0);", """      if constexpr (requires { model.setExternalTime(tMinust0); })
        model.setExternalTime(tMinust0);

""")
path = output / relative
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(text)
