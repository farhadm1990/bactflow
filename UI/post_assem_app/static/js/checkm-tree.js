(function (global) {
  "use strict";

  const RANK_COLORS = {
    d: "#1f4e5f",
    k: "#1f4e5f",
    p: "#2a6f7f",
    c: "#356f62",
    o: "#42756b",
    f: "#548a68",
    g: "#6aa37a",
    s: "#3d5c9a"
  };

  const state = {
    newick: "",
    layout: "rectangular",
    tree: null,
    svg: null
  };

  function parseTaxon(name) {
    const raw = String(name || "").trim();
    const m = raw.match(/^([a-z])__(.+)$/i);
    if (m) {
      return {
        rank: m[1].toLowerCase(),
        label: m[2].replace(/_/g, " ").trim()
      };
    }
    return { rank: "", label: raw.replace(/_/g, " ").trim() };
  }

  function parseNewick(text) {
    const s = String(text || "").replace(/\[[^\]]*\]/g, "").trim().replace(/;$/, "");
    if (!s) {
      throw new Error("empty newick");
    }
    const root = { name: "", length: null, children: [] };
    let current = root;
    const stack = [];
    let buf = "";
    let inQuote = false;
    let quoteChar = "";

    function flush(node) {
      const chunk = buf.trim();
      buf = "";
      if (!chunk) {
        return;
      }
      let name = chunk;
      let length = null;
      if (name[0] === "'" || name[0] === '"') {
        const q = name[0];
        let j = 1;
        while (j < name.length) {
          if (name[j] === q) {
            if (name[j + 1] === q) {
              j += 2;
              continue;
            }
            const inner = name.slice(1, j).replace(new RegExp(q + q, "g"), q);
            const rest = name.slice(j + 1);
            name = inner;
            if (rest.charAt(0) === ":") {
              const n = parseFloat(rest.slice(1));
              length = Number.isFinite(n) ? n : null;
            }
            break;
          }
          j += 1;
        }
      } else if (name.includes(":")) {
        const idx = name.lastIndexOf(":");
        const dist = parseFloat(name.slice(idx + 1));
        name = name.slice(0, idx);
        length = Number.isFinite(dist) ? dist : null;
      }
      node.name = name;
      node.length = length;
    }

    for (let i = 0; i < s.length; i += 1) {
      const c = s[i];
      if (inQuote) {
        buf += c;
        if (c === quoteChar) {
          inQuote = false;
        }
        continue;
      }
      if (c === "'" || c === '"') {
        inQuote = true;
        quoteChar = c;
        buf += c;
        continue;
      }
      if (c === "(") {
        const child = { name: "", length: null, children: [] };
        current.children.push(child);
        stack.push(current);
        current = child;
      } else if (c === ",") {
        flush(current);
        const parent = stack[stack.length - 1];
        const child = { name: "", length: null, children: [] };
        parent.children.push(child);
        current = child;
      } else if (c === ")") {
        flush(current);
        current = stack.pop();
      } else {
        buf += c;
      }
    }
    flush(current);
    return root;
  }

  function walk(node, fn) {
    fn(node);
    (node.children || []).forEach((child) => walk(child, fn));
  }

  function tipsOf(node) {
    const tips = [];
    walk(node, (n) => {
      if (!n.children || !n.children.length) {
        tips.push(n);
      }
    });
    return tips;
  }

  function maxDepthDist(node, acc) {
    const here = acc + (Number(node.length) || 0);
    if (!node.children || !node.children.length) {
      return here;
    }
    return Math.max.apply(null, node.children.map((c) => maxDepthDist(c, here)));
  }

  function formatLen(value) {
    const n = Number(value);
    if (!Number.isFinite(n)) {
      return "";
    }
    const abs = Math.abs(n);
    if (abs === 0) {
      return "0";
    }
    if (abs >= 1) {
      return n.toFixed(3).replace(/\.?0+$/, "");
    }
    if (abs >= 0.001) {
      return n.toFixed(4).replace(/\.?0+$/, "");
    }
    return n.toExponential(2);
  }

  function esc(value) {
    return String(value ?? "")
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  function layoutRectangular(root, width) {
    const tips = tipsOf(root);
    const row = 34;
    const pad = { l: 36, r: 28, t: 48, b: 36 };
    const longest = tips.reduce((m, t) => Math.max(m, (parseTaxon(t.name).label || "").length), 8);
    pad.r = Math.max(140, Math.min(320, longest * 8 + 24));
    const height = Math.max(420, pad.t + pad.b + tips.length * row);
    const innerW = Math.max(220, width - pad.l - pad.r);
    const innerH = height - pad.t - pad.b;
    const maxD = maxDepthDist(root, 0) || 1;
    tips.forEach((t, i) => {
      t.y = pad.t + ((i + 0.5) * innerH) / tips.length;
    });
    function setY(n) {
      if (n.children && n.children.length) {
        n.children.forEach(setY);
        n.y = (n.children[0].y + n.children[n.children.length - 1].y) / 2;
      }
    }
    setY(root);
    function setX(n, acc) {
      const here = acc + (Number(n.length) || 0);
      n.x = pad.l + (here / maxD) * innerW;
      (n.children || []).forEach((c) => setX(c, here));
    }
    setX(root, 0);
    return { width, height, pad, maxD };
  }

  function layoutCircular(root, size) {
    const tips = tipsOf(root);
    const pad = 28;
    const longest = tips.reduce((m, t) => Math.max(m, (parseTaxon(t.name).label || "").length), 8);
    const labelR = Math.max(70, Math.min(180, longest * 6.2));
    const cx = size / 2;
    const cy = size / 2;
    const maxR = Math.max(80, size / 2 - pad - labelR);
    const maxD = maxDepthDist(root, 0) || 1;
    tips.forEach((t, i) => {
      t.angle = -Math.PI / 2 + (i / Math.max(tips.length, 1)) * Math.PI * 2;
      if (tips.length === 1) {
        t.angle = -Math.PI / 2;
      }
    });
    function setAngle(n) {
      if (n.children && n.children.length) {
        n.children.forEach(setAngle);
        n.angle = (n.children[0].angle + n.children[n.children.length - 1].angle) / 2;
      }
    }
    setAngle(root);
    function setR(n, acc) {
      const here = acc + (Number(n.length) || 0);
      n.r = 16 + (here / maxD) * (maxR - 16);
      n.x = cx + n.r * Math.cos(n.angle);
      n.y = cy + n.r * Math.sin(n.angle);
      (n.children || []).forEach((c) => setR(c, here));
    }
    setR(root, 0);
    root.x = cx;
    root.y = cy;
    return { width: size, height: size, cx, cy, maxR, labelR, maxD };
  }

  function taxonBox(x, y, name, opts) {
    const parsed = parseTaxon(name);
    if (!parsed.label) {
      return "";
    }
    const rank = parsed.rank;
    const color = RANK_COLORS[rank] || "#42756b";
    const text = rank ? `${rank}  ${parsed.label}` : parsed.label;
    const font = 11;
    const w = Math.min(240, 18 + text.length * 6.35);
    const h = 18;
    const rx = x - w / 2;
    const ry = y - h / 2 + (opts && opts.dy ? opts.dy : 0);
    return `<g class="checkm-taxon">
      <rect x="${rx.toFixed(1)}" y="${ry.toFixed(1)}" width="${w.toFixed(1)}" height="${h}" rx="4" ry="4"
        fill="${color}" stroke="#1b2e2a" stroke-width="0.6"/>
      <text x="${x.toFixed(1)}" y="${(ry + 13).toFixed(1)}" text-anchor="middle"
        font-family="Georgia, 'Times New Roman', serif" font-size="${font}" fill="#fff">${esc(text)}</text>
    </g>`;
  }

  function drawRectangular(root, meta) {
    const parts = [];
    parts.push(`<rect width="${meta.width}" height="${meta.height}" fill="#f7faf8"/>`);
    walk(root, (n) => {
      if (!n.children || !n.children.length) {
        return;
      }
      const y0 = n.children[0].y;
      const y1 = n.children[n.children.length - 1].y;
      parts.push(`<line x1="${n.x.toFixed(1)}" y1="${y0.toFixed(1)}" x2="${n.x.toFixed(1)}" y2="${y1.toFixed(1)}" stroke="#2f4f4f" stroke-width="1.6"/>`);
      n.children.forEach((c) => {
        parts.push(`<line x1="${n.x.toFixed(1)}" y1="${c.y.toFixed(1)}" x2="${c.x.toFixed(1)}" y2="${c.y.toFixed(1)}" stroke="#2f4f4f" stroke-width="1.6"/>`);
        const dx = c.x - n.x;
        if (dx > 28 && c.length != null) {
          const mx = n.x + dx / 2;
          parts.push(`<text x="${mx.toFixed(1)}" y="${(c.y - 4).toFixed(1)}" text-anchor="middle"
            font-family="ui-monospace, Menlo, monospace" font-size="9" fill="#6a8078">${esc(formatLen(c.length))}</text>`);
        }
      });
    });
    walk(root, (n) => {
      if (n.children && n.children.length) {
        if (n.name) {
          parts.push(taxonBox(n.x, n.y, n.name, { dy: -14 }));
        } else {
          parts.push(`<circle cx="${n.x.toFixed(1)}" cy="${n.y.toFixed(1)}" r="3.2" fill="#42756b" stroke="#fff" stroke-width="1"/>`);
        }
        return;
      }
      const tip = parseTaxon(n.name);
      parts.push(`<circle cx="${n.x.toFixed(1)}" cy="${n.y.toFixed(1)}" r="4" fill="#c45c26" stroke="#fff" stroke-width="1"/>`);
      parts.push(`<text x="${(n.x + 10).toFixed(1)}" y="${(n.y + 4).toFixed(1)}"
        font-family="Georgia, 'Times New Roman', serif" font-size="13" font-style="italic" fill="#1b2e2a">${esc(tip.label)}</text>`);
    });
    parts.push(`<text x="16" y="22" font-family="system-ui, sans-serif" font-size="11" fill="#5b6e68">Rectangular tree · boxed names are taxon ranks · italic tips are species · numbers are branch lengths</text>`);
    return parts.join("");
  }

  function polarLine(x1, y1, x2, y2) {
    return `<line x1="${x1.toFixed(1)}" y1="${y1.toFixed(1)}" x2="${x2.toFixed(1)}" y2="${y2.toFixed(1)}" stroke="#2f4f4f" stroke-width="1.5"/>`;
  }

  function arcPath(cx, cy, r, a0, a1) {
    if (Math.abs(a1 - a0) < 0.002 || r < 1) {
      return "";
    }
    const x0 = cx + r * Math.cos(a0);
    const y0 = cy + r * Math.sin(a0);
    const x1 = cx + r * Math.cos(a1);
    const y1 = cy + r * Math.sin(a1);
    const large = Math.abs(a1 - a0) > Math.PI ? 1 : 0;
    const sweep = a1 > a0 ? 1 : 0;
    return `<path d="M ${x0.toFixed(1)} ${y0.toFixed(1)} A ${r.toFixed(1)} ${r.toFixed(1)} 0 ${large} ${sweep} ${x1.toFixed(1)} ${y1.toFixed(1)}" fill="none" stroke="#2f4f4f" stroke-width="1.5"/>`;
  }

  function drawCircular(root, meta) {
    const parts = [];
    parts.push(`<rect width="${meta.width}" height="${meta.height}" fill="#f7faf8"/>`);
    walk(root, (n) => {
      if (!n.children || !n.children.length) {
        return;
      }
      const a0 = n.children[0].angle;
      const a1 = n.children[n.children.length - 1].angle;
      parts.push(arcPath(meta.cx, meta.cy, n.r, a0, a1));
      n.children.forEach((c) => {
        const px = meta.cx + n.r * Math.cos(c.angle);
        const py = meta.cy + n.r * Math.sin(c.angle);
        parts.push(polarLine(px, py, c.x, c.y));
        const mx = (px + c.x) / 2;
        const my = (py + c.y) / 2;
        const dist = Math.hypot(c.x - px, c.y - py);
        if (dist > 26 && c.length != null) {
          parts.push(`<text x="${mx.toFixed(1)}" y="${my.toFixed(1)}" text-anchor="middle"
            font-family="ui-monospace, Menlo, monospace" font-size="9" fill="#6a8078">${esc(formatLen(c.length))}</text>`);
        }
      });
    });
    walk(root, (n) => {
      if (n.children && n.children.length) {
        if (n.name) {
          parts.push(taxonBox(n.x, n.y, n.name, {}));
        } else {
          parts.push(`<circle cx="${n.x.toFixed(1)}" cy="${n.y.toFixed(1)}" r="3" fill="#42756b" stroke="#fff" stroke-width="1"/>`);
        }
        return;
      }
      const tip = parseTaxon(n.name);
      const angle = n.angle;
      let deg = (angle * 180) / Math.PI;
      let anchor = "start";
      const lx = n.x + Math.cos(angle) * 10;
      const ly = n.y + Math.sin(angle) * 10;
      if (angle > Math.PI / 2 || angle < -Math.PI / 2) {
        deg += 180;
        anchor = "end";
      }
      parts.push(`<circle cx="${n.x.toFixed(1)}" cy="${n.y.toFixed(1)}" r="4" fill="#c45c26" stroke="#fff" stroke-width="1"/>`);
      parts.push(`<text x="${lx.toFixed(1)}" y="${ly.toFixed(1)}" text-anchor="${anchor}" dominant-baseline="middle"
        transform="rotate(${deg.toFixed(2)} ${lx.toFixed(1)} ${ly.toFixed(1)})"
        font-family="Georgia, 'Times New Roman', serif" font-size="12" font-style="italic" fill="#1b2e2a">${esc(tip.label)}</text>`);
    });
    parts.push(`<text x="16" y="22" font-family="system-ui, sans-serif" font-size="11" fill="#5b6e68">Circular tree · boxed names are taxon ranks · italic tips are species · numbers are branch lengths</text>`);
    return parts.join("");
  }

  function renderInto(svg, newick, layout) {
    state.newick = newick;
    state.layout = layout || state.layout || "rectangular";
    state.svg = svg;
    const wrap = svg.parentElement;
    const width = Math.max(640, (wrap && wrap.clientWidth) || 900);
    let tree;
    try {
      tree = parseNewick(newick);
    } catch (_err) {
      svg.innerHTML = "";
      svg.setAttribute("viewBox", "0 0 640 120");
      svg.innerHTML = `<rect width="640" height="120" fill="#f7faf8"/><text x="24" y="64" fill="#a33">Could not parse taxon_tree.newick</text>`;
      return;
    }
    state.tree = tree;
    let inner = "";
    let w = width;
    let h = 480;
    if (state.layout === "circular") {
      const size = Math.max(720, Math.min(1200, width - 8));
      const meta = layoutCircular(tree, size);
      inner = drawCircular(tree, meta);
      w = meta.width;
      h = meta.height;
    } else {
      const meta = layoutRectangular(tree, width);
      inner = drawRectangular(tree, meta);
      w = meta.width;
      h = meta.height;
    }
    svg.setAttribute("viewBox", `0 0 ${w} ${h}`);
    svg.setAttribute("width", String(w));
    svg.setAttribute("height", String(h));
    svg.innerHTML = inner;
  }

  function setLayout(layout) {
    if (!state.svg || !state.newick) {
      return;
    }
    renderInto(state.svg, state.newick, layout);
  }

  function downloadPng(filename) {
    const svg = state.svg;
    if (!svg) {
      return;
    }
    const w = svg.viewBox.baseVal.width || svg.clientWidth || 900;
    const h = svg.viewBox.baseVal.height || svg.clientHeight || 640;
    const clone = svg.cloneNode(true);
    clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    const xml = new XMLSerializer().serializeToString(clone);
    const blob = new Blob([xml], { type: "image/svg+xml;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const img = new Image();
    img.onload = () => {
      const scale = 2;
      const canvas = document.createElement("canvas");
      canvas.width = Math.max(1, Math.round(w * scale));
      canvas.height = Math.max(1, Math.round(h * scale));
      const ctx = canvas.getContext("2d");
      ctx.fillStyle = "#f7faf8";
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
      URL.revokeObjectURL(url);
      canvas.toBlob((png) => {
        if (!png) {
          return;
        }
        const a = document.createElement("a");
        a.href = URL.createObjectURL(png);
        a.download = filename || "checkm_taxon_tree.png";
        document.body.appendChild(a);
        a.click();
        a.remove();
      }, "image/png");
    };
    img.onerror = () => URL.revokeObjectURL(url);
    img.src = url;
  }

  global.CheckmTreeViz = {
    render: renderInto,
    setLayout,
    downloadPng,
    getLayout: () => state.layout
  };
})(window);
