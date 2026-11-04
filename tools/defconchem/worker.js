/* defconchem web worker: runs the Python engine (py/defconchem.py) in Pyodide,
   so the page stays responsive and no data ever leaves the user's computer. */
const PYODIDE_VERSION = "0.27.2";
// Default: public CDN. Self-hosting (e.g. offline lab server): open the page with
// ?pyodide=<URL of a folder containing pyodide.js, its wasm/stdlib and the numpy/scipy/openblas packages>
const PYODIDE_BASE = (() => {
  const b = new URL(self.location.href).searchParams.get("base");
  return b ? (b.endsWith("/") ? b : b + "/") : `https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/full/`;
})();
importScripts(PYODIDE_BASE + "pyodide.js");

let pyodide = null;

self.reportProgress = (done, total) => postMessage({ type: "progress", done, total });

async function init() {
  try {
    postMessage({ type: "status", text: "Loading Python runtime…" });
    pyodide = await loadPyodide({ indexURL: PYODIDE_BASE });
    postMessage({ type: "status", text: "Loading numpy and scipy…" });
    await pyodide.loadPackage(["numpy", "scipy"]);
    postMessage({ type: "status", text: "Loading defconchem…" });
    const src = await (await fetch("py/defconchem.py", { cache: "no-cache" })).text();
    pyodide.FS.mkdirTree("/work");
    pyodide.FS.writeFile("/home/pyodide/defconchem.py", src);
    await pyodide.runPythonAsync(`
import sys, io, json, math, os, contextlib, shutil
sys.path.insert(0, "/home/pyodide")
import defconchem
import js

def _clean(x):
    if isinstance(x, float):
        return x if math.isfinite(x) else None
    if isinstance(x, dict):
        return {k: _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    return x

def web_run(check):
    os.chdir("/work")
    out, err = io.StringIO(), io.StringIO()
    res, ok = None, True
    try:
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            res = defconchem.run("input.txt", check=check, plot=False, write_files=True,
                                 progress=lambda d, t: js.reportProgress(d, t))
    except SystemExit:
        ok = False
    except Exception as e:
        ok = False
        err.write("\\nERROR: " + type(e).__name__ + ": " + str(e) + "\\n")
    files = {}
    if res:
        for f in res.get("files", []):
            if os.path.exists(f):
                files[f] = open(f).read()
    return json.dumps({"ok": ok, "log": out.getvalue(), "err": err.getvalue(),
                       "result": _clean(res), "files": files})

def web_template():
    return defconchem.TEMPLATE

def web_version():
    return defconchem.__version__
`);
    const version = pyodide.runPython("web_version()");
    const template = pyodide.runPython("web_template()");
    postMessage({ type: "ready", version, template });
  } catch (e) {
    postMessage({ type: "fatal", text: String(e) });
  }
}

const ready = init();

self.onmessage = async (ev) => {
  await ready;
  const { cmd, input, files, check, id } = ev.data;
  if (cmd !== "run") return;
  try {
    // fresh working directory for every run
    pyodide.runPython(`
import os, shutil
os.chdir("/")
shutil.rmtree("/work", ignore_errors=True)
os.makedirs("/work", exist_ok=True)
for _f in os.listdir("/work"):
    _p = os.path.join("/work", _f)
    shutil.rmtree(_p, ignore_errors=True) if os.path.isdir(_p) else os.remove(_p)`);
    for (const f of files || []) {
      pyodide.FS.writeFile(`/work/${f.name}`, new Uint8Array(f.data));
    }
    pyodide.FS.writeFile("/work/input.txt", input);
    const t0 = performance.now();
    const json = pyodide.runPython(`web_run(${check ? "True" : "False"})`);
    const msg = JSON.parse(json);
    msg.seconds = (performance.now() - t0) / 1000;
    postMessage({ type: "result", id, ...msg });
  } catch (e) {
    postMessage({ type: "result", id, ok: false, log: "", err: String(e), result: null, files: {} });
  }
};
