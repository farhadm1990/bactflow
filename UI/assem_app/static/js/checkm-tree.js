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

  const instances = {};
  let lastSvg = null;

  function instanceOf(svg) {
    const node = svg || lastSvg;
    const key = (node && node.id) || "default";
    if (!instances[key]) {
      instances[key] = { newick: "", layout: "rectangular", tree: null, svg: null };
    }
    return instances[key];
  }

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

  let measureCtx = null;
  function measureTextWidth(text, font) {
    if (typeof document !== "undefined") {
      if (!measureCtx) {
        measureCtx = document.createElement("canvas").getContext("2d");
      }
      if (measureCtx) {
        measureCtx.font = font;
        return measureCtx.measureText(text).width;
      }
    }
    return String(text || "").length * 7.1;
  }

  function isNamedInternal(n) {
    return !!(n && n.children && n.children.length && parseTaxon(n.name).label);
  }

  function maxNamedChain(n) {
    const here = isNamedInternal(n) ? 1 : 0;
    if (!n.children || !n.children.length) {
      return here;
    }
    return here + Math.max.apply(null, n.children.map(maxNamedChain));
  }

  function taxonStyle(name) {
    const parsed = parseTaxon(name);
    const rank = parsed.rank;
    const text = rank && parsed.label ? `${rank}  ${parsed.label}` : parsed.label;
    const font = "12px Georgia, 'Times New Roman', serif";
    const w = Math.ceil(Math.min(320, Math.max(52, measureTextWidth(text, font) + 18)));
    const h = 22;
    return {
      parsed,
      rank,
      text,
      color: RANK_COLORS[rank] || "#42756b",
      w,
      h,
      font
    };
  }

  function boxesOverlap(a, b, gap) {
    const g = gap == null ? 8 : gap;
    return Math.abs(a.x - b.x) < (a.w + b.w) / 2 + g && Math.abs(a.y - b.y) < (a.h + b.h) / 2 + g;
  }

  function resolveLabelCollisions(boxes, bounds) {
    const pad = 8;
    for (let iter = 0; iter < 90; iter += 1) {
      let moved = false;
      for (let i = 0; i < boxes.length; i += 1) {
        for (let j = i + 1; j < boxes.length; j += 1) {
          const a = boxes[i];
          const b = boxes[j];
          if (!boxesOverlap(a, b, 10)) {
            continue;
          }
          let dx = a.x - b.x;
          let dy = a.y - b.y;
          if (dx === 0 && dy === 0) {
            dx = 0.6;
            dy = (i % 2 === 0 ? 1 : -1) * 0.8;
          }
          const len = Math.hypot(dx, dy) || 1;
          const overlapX = (a.w + b.w) / 2 + 10 - Math.abs(a.x - b.x);
          const overlapY = (a.h + b.h) / 2 + 10 - Math.abs(a.y - b.y);
          const push = Math.max(1.5, Math.min(overlapX, overlapY) / 2 + 1);
          const ux = dx / len;
          const uy = dy / len;
          const wa = a.kind === "tip" ? 0.12 : 1;
          const wb = b.kind === "tip" ? 0.12 : 1;
          a.x += ux * push * wa;
          a.y += uy * push * wa;
          b.x -= ux * push * wb;
          b.y -= uy * push * wb;
          moved = true;
        }
      }
      boxes.forEach((box) => {
        const hw = box.w / 2;
        const hh = box.h / 2;
        const nx = Math.min(bounds.x1 - hw - pad, Math.max(bounds.x0 + hw + pad, box.x));
        const ny = Math.min(bounds.y1 - hh - pad, Math.max(bounds.y0 + hh + pad, box.y));
        if (nx !== box.x || ny !== box.y) {
          box.x = nx;
          box.y = ny;
          moved = true;
        }
      });
      if (!moved) {
        break;
      }
    }
  }

  function layoutRectangular(root, width) {
    const tips = tipsOf(root);
    const named = maxNamedChain(root);
    const row = Math.max(40, 28 + Math.min(24, named));
    const pad = { l: 28, r: 28, t: 52, b: 40 };
    const longest = tips.reduce((m, t) => Math.max(m, (parseTaxon(t.name).label || "").length), 8);
    pad.r = Math.max(160, Math.min(360, longest * 8.2 + 28));
    pad.l = 36;
    const minInnerW = 180 + named * 96;
    const height = Math.max(460, pad.t + pad.b + Math.max(tips.length, 2) * row);
    const innerW = Math.max(minInnerW, width - pad.l - pad.r);
    const innerH = height - pad.t - pad.b;
    const maxD = maxDepthDist(root, 0) || 1;
    const minNamed = 88;
    const minChild = 36;
    tips.forEach((t, i) => {
      t.y = pad.t + ((i + 0.5) * innerH) / Math.max(tips.length, 1);
    });
    function setY(n) {
      if (n.children && n.children.length) {
        n.children.forEach(setY);
        n.y = (n.children[0].y + n.children[n.children.length - 1].y) / 2;
      }
    }
    setY(root);
    function setX(n, parentX) {
      if (n === root) {
        n.x = pad.l;
      } else {
        const natural = ((Number(n.length) || 0) / maxD) * innerW;
        const minStep = isNamedInternal(n) ? minNamed : minChild;
        n.x = parentX + Math.max(natural, minStep);
      }
      (n.children || []).forEach((c) => setX(c, n.x));
    }
    setX(root, pad.l);
    let maxX = pad.l;
    walk(root, (n) => {
      maxX = Math.max(maxX, n.x);
    });
    const totalW = pad.l + Math.max(innerW, maxX - pad.l + 8) + pad.r;
    return { width: totalW, height, pad, maxD };
  }

  function layoutCircular(root, size) {
    const tips = tipsOf(root);
    const named = maxNamedChain(root);
    const longest = tips.reduce((m, t) => Math.max(m, (parseTaxon(t.name).label || "").length), 8);
    const labelR = Math.max(90, Math.min(220, longest * 7.2));
    const minNamed = 56;
    const minChild = 30;
    const needR = 24 + named * minNamed + Math.max(tips.length, 1) * 8 + 70;
    size = Math.max(size, Math.ceil((needR + labelR + 80) * 2));
    const cx = size / 2;
    const cy = size / 2;
    const maxR = Math.max(120, size / 2 - 24 - labelR);
    const maxD = maxDepthDist(root, 0) || 1;
    const nTips = Math.max(tips.length, 1);
    const span = Math.PI * 2;
    const start = -Math.PI / 2;
    tips.forEach((t, i) => {
      t.angle = nTips === 1 ? -Math.PI / 2 : start + (i / nTips) * span;
    });
    function setAngle(n) {
      if (n.children && n.children.length) {
        n.children.forEach(setAngle);
        n.angle = (n.children[0].angle + n.children[n.children.length - 1].angle) / 2;
      }
    }
    setAngle(root);
    function setR(n, parentR) {
      if (n === root) {
        n.r = 0;
      } else {
        const natural = ((Number(n.length) || 0) / maxD) * maxR;
        const minStep = isNamedInternal(n) ? minNamed : minChild;
        n.r = parentR + Math.max(natural, minStep);
      }
      n.x = cx + n.r * Math.cos(n.angle || 0);
      n.y = cy + n.r * Math.sin(n.angle || 0);
      (n.children || []).forEach((c) => setR(c, n.r));
    }
    root.angle = root.angle || -Math.PI / 2;
    setR(root, 0);
    root.x = cx;
    root.y = cy;
    return { width: size, height: size, cx, cy, maxR, labelR, maxD };
  }

  function drawTaxonBox(box) {
    const rx = box.x - box.w / 2;
    const ry = box.y - box.h / 2;
    let leader = "";
    const dist = Math.hypot(box.x - box.ax, box.y - box.ay);
    if (dist > 14) {
      leader = `<line x1="${box.ax.toFixed(1)}" y1="${box.ay.toFixed(1)}" x2="${box.x.toFixed(1)}" y2="${box.y.toFixed(1)}"
        stroke="#42756b" stroke-width="0.9" stroke-dasharray="3 2"/>`;
    }
    return `${leader}<g class="checkm-taxon">
      <rect x="${rx.toFixed(1)}" y="${ry.toFixed(1)}" width="${box.w}" height="${box.h}" rx="5" ry="5"
        fill="${box.color}" stroke="#1b2e2a" stroke-width="0.7"/>
      <text x="${box.x.toFixed(1)}" y="${(box.y + 4.5).toFixed(1)}" text-anchor="middle"
        font-family="Georgia, 'Times New Roman', serif" font-size="12" fill="#fff">${esc(box.text)}</text>
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
        if (dx > 36 && c.length != null) {
          const mx = n.x + dx / 2;
          parts.push(`<text x="${mx.toFixed(1)}" y="${(c.y - 5).toFixed(1)}" text-anchor="middle"
            font-family="ui-monospace, Menlo, monospace" font-size="9" fill="#6a8078">${esc(formatLen(c.length))}</text>`);
        }
      });
    });
    const labels = [];
    walk(root, (n) => {
      if (n.children && n.children.length) {
        if (!parseTaxon(n.name).label) {
          parts.push(`<circle cx="${n.x.toFixed(1)}" cy="${n.y.toFixed(1)}" r="3.2" fill="#42756b" stroke="#fff" stroke-width="1"/>`);
          return;
        }
        const st = taxonStyle(n.name);
        labels.push({
          kind: "taxon",
          ax: n.x,
          ay: n.y,
          x: n.x,
          y: n.y - 20,
          w: st.w,
          h: st.h,
          text: st.text,
          color: st.color
        });
        return;
      }
      const tip = parseTaxon(n.name);
      const st = taxonStyle(tip.label);
      labels.push({
        kind: "tip",
        ax: n.x,
        ay: n.y,
        x: n.x + 12 + st.w / 2,
        y: n.y,
        w: st.w,
        h: 18,
        text: tip.label,
        color: null,
        tip: true
      });
      parts.push(`<circle cx="${n.x.toFixed(1)}" cy="${n.y.toFixed(1)}" r="4" fill="#c45c26" stroke="#fff" stroke-width="1"/>`);
    });
    resolveLabelCollisions(labels, { x0: 4, y0: 28, x1: meta.width - 4, y1: meta.height - 4 });
    labels.forEach((box) => {
      if (box.tip) {
        parts.push(`<text x="${(box.x - box.w / 2).toFixed(1)}" y="${(box.y + 4).toFixed(1)}"
          font-family="Georgia, 'Times New Roman', serif" font-size="13" font-style="italic" fill="#1b2e2a">${esc(box.text)}</text>`);
        return;
      }
      parts.push(drawTaxonBox(box));
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
        if (dist > 34 && c.length != null) {
          parts.push(`<text x="${mx.toFixed(1)}" y="${my.toFixed(1)}" text-anchor="middle"
            font-family="ui-monospace, Menlo, monospace" font-size="9" fill="#6a8078">${esc(formatLen(c.length))}</text>`);
        }
      });
    });
    const labels = [];
    let side = 1;
    walk(root, (n) => {
      if (n.children && n.children.length) {
        parts.push(`<circle cx="${n.x.toFixed(1)}" cy="${n.y.toFixed(1)}" r="3.2" fill="#42756b" stroke="#fff" stroke-width="1"/>`);
        if (!parseTaxon(n.name).label) {
          return;
        }
        const st = taxonStyle(n.name);
        const ang = n.angle || 0;
        const perp = ang + Math.PI / 2;
        side *= -1;
        const along = n === root ? 34 : 10;
        const beside = n === root ? 0 : side * (st.w / 2 + 16);
        labels.push({
          kind: "taxon",
          ax: n.x,
          ay: n.y,
          x: n.x + Math.cos(ang) * along + Math.cos(perp) * beside,
          y: n.y + Math.sin(ang) * along + Math.sin(perp) * beside,
          w: st.w,
          h: st.h,
          text: st.text,
          color: st.color
        });
        return;
      }
      const tip = parseTaxon(n.name);
      const angle = n.angle;
      let deg = (angle * 180) / Math.PI;
      let anchor = "start";
      const lx = n.x + Math.cos(angle) * 12;
      const ly = n.y + Math.sin(angle) * 12;
      if (angle > Math.PI / 2 || angle < -Math.PI / 2) {
        deg += 180;
        anchor = "end";
      }
      parts.push(`<circle cx="${n.x.toFixed(1)}" cy="${n.y.toFixed(1)}" r="4" fill="#c45c26" stroke="#fff" stroke-width="1"/>`);
      parts.push(`<text x="${lx.toFixed(1)}" y="${ly.toFixed(1)}" text-anchor="${anchor}" dominant-baseline="middle"
        transform="rotate(${deg.toFixed(2)} ${lx.toFixed(1)} ${ly.toFixed(1)})"
        font-family="Georgia, 'Times New Roman', serif" font-size="12" font-style="italic" fill="#1b2e2a">${esc(tip.label)}</text>`);
    });
    resolveLabelCollisions(labels, { x0: 8, y0: 28, x1: meta.width - 8, y1: meta.height - 8 });
    labels.forEach((box) => parts.push(drawTaxonBox(box)));
    parts.push(`<text x="16" y="22" font-family="system-ui, sans-serif" font-size="11" fill="#5b6e68">Circular tree · boxed names are taxon ranks · italic tips are species · numbers are branch lengths</text>`);
    return parts.join("");
  }

  function renderInto(svg, newick, layout) {
    lastSvg = svg;
    const inst = instanceOf(svg);
    inst.newick = newick;
    inst.layout = layout || inst.layout || "rectangular";
    inst.svg = svg;
    const wrap = svg.parentElement;
    const width = Math.max(640, (wrap && wrap.clientWidth) || 900);
    let tree;
    try {
      tree = parseNewick(newick);
    } catch (_err) {
      svg.innerHTML = "";
      svg.setAttribute("viewBox", "0 0 640 120");
      svg.innerHTML = `<rect width="640" height="120" fill="#f7faf8"/><text x="24" y="64" fill="#a33">Could not parse tree</text>`;
      return;
    }
    inst.tree = tree;
    let inner = "";
    let w = width;
    let h = 480;
    if (inst.layout === "circular") {
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

  function setLayout(layout, svg) {
    const inst = instanceOf(svg);
    if (!inst.svg || !inst.newick) {
      return;
    }
    renderInto(inst.svg, inst.newick, layout);
  }

  function downloadPng(filename, svg) {
    const inst = instanceOf(svg);
    const target = inst.svg;
    if (!target) {
      return;
    }
    const w = target.viewBox.baseVal.width || target.clientWidth || 900;
    const h = target.viewBox.baseVal.height || target.clientHeight || 640;
    const clone = target.cloneNode(true);
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
    getLayout: (svg) => instanceOf(svg).layout
  };
})(window);
