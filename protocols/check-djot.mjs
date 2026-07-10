import { parse, renderHTML } from "@djot/djot";
import { readFileSync } from "fs";
let bad = 0;
for (const f of process.argv.slice(2)) {
  try {
    let warned = false;
    parse(readFileSync(f, "utf8"), { warn: (w) => { warned = true; console.log("WARN", f, w.message ?? String(w)); } });
    renderHTML(parse(readFileSync(f, "utf8")));
    if (!warned) console.log("OK  ", f);
  } catch (e) { bad++; console.log("FAIL", f, e.message); }
}
process.exit(bad ? 1 : 0);
