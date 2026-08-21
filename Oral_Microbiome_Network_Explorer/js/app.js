/* Oral Microbiome Network Explorer
   Runs entirely in the browser and is suitable for GitHub Pages.
*/

const state = {
  currentSubsite: APP_CONFIG.defaultSubsite,
  activeSubsites: [APP_CONFIG.defaultSubsite],
  selectedSpecies: null,
  crossMode: false,
  crossClusterFilters: {},
  cache: {},
  cy: null
};

const $ = (id) => document.getElementById(id);

function speciesFromRow(row) {
  // The current PCAdf CSVs store species names in an unnamed first column.
  // PapaParse exposes that column with the empty-string key "".
  // Keep fallbacks so the app also works if the files are later saved with
  // an explicit species/taxon column name.
  return String(
    row[""] ??
    row["Unnamed: 0"] ??
    row.species ??
    row.Species ??
    row.taxon ??
    row.Taxon ??
    ""
  ).trim();
}

function cleanSpeciesName(name) {
  return String(name).replaceAll("_", " ");
}

function fmt(v) {
  if (v === null || v === undefined || v === "") return "—";
  if (typeof v === "number") {
    if (!Number.isFinite(v)) return "—";
    return Math.abs(v) < 100 ? v.toFixed(3).replace(/0+$/,"").replace(/\.$/,"") : v.toLocaleString();
  }
  return String(v);
}

function titleCaseField(s) {
  const exact = {
    cluster: "Cluster",
    CS: "CS",
    HS: "HS",
    SS: "SS",
    sHACK: "sHACK",
    HAC: "HAC"
  };
  if (exact[s]) return exact[s];
  return String(s)
    .replaceAll("_", " ")
    .replace(/\b\w/g, m => m.toUpperCase());
}

function isHiddenMetadataField(field) {
  return APP_CONFIG.hiddenMetadataPatterns.some(re => re.test(field));
}

function papa(url) {
  return new Promise((resolve, reject) => {
    Papa.parse(url, {
      download: true,
      header: true,
      dynamicTyping: true,
      skipEmptyLines: true,
      complete: (r) => resolve(r.data),
      error: reject
    });
  });
}

async function loadSubsite(key) {
  if (state.cache[key]) return state.cache[key];

  const cfg = APP_CONFIG.subsites[key];
  const [nodeRows, edgeRows] = await Promise.all([papa(cfg.nodeFile), papa(cfg.edgeFile)]);

  const nodeMap = new Map();
  nodeRows.forEach(row => {
    const species = speciesFromRow(row);
    if (species) nodeMap.set(species, { ...row, species });
  });

  const edges = edgeRows
    .map((row, i) => ({
      id: `${key}::e${i}`,
      source: String(row.Node1 ?? row.source ?? "").trim(),
      target: String(row.Node2 ?? row.target ?? "").trim(),
      weight: Number(row.Weight ?? row.weight ?? 1)
    }))
    .filter(e => e.source && e.target);

  const connected = new Set();
  edges.forEach(e => { connected.add(e.source); connected.add(e.target); });

  const data = { nodeRows, nodeMap, edges, connected };
  state.cache[key] = data;
  return data;
}

async function preloadAll() {
  await Promise.all(Object.keys(APP_CONFIG.subsites).map(loadSubsite));
}

function clusterColor(subsite, cluster) {
  return APP_CONFIG.subsites[subsite].clusterColors[String(cluster)] || "#94a3b8";
}

function mainScoreInfo(subsite, row) {
  const field = subsite === "saliva" ? "sHACK" : "HAC";
  const value = Number(row?.[field]);
  const high = Number.isFinite(value) && value >= 0.90;
  return {
    field,
    value,
    high,
    label: subsite === "saliva" ? "sHACK" : "HAC",
    highLabel: subsite === "saliva" ? "High sHACK" : "High HAC"
  };
}

function neighborRowHtml(subsite, species) {
  const row = state.cache[subsite]?.nodeMap.get(species) || {};
  const score = mainScoreInfo(subsite, row);
  const scoreText = Number.isFinite(score.value) ? score.value.toFixed(3) : "—";
  const highBadge = score.high ? `<span class="high-score-badge">${score.highLabel}</span>` : "";
  return `<button class="neighbor" data-species="${escapeHtml(species)}">
    <span class="neighbor-name">${escapeHtml(cleanSpeciesName(species))}</span>
    <span class="neighbor-details">
      <span class="neighbor-cluster">C${escapeHtml(String(row.cluster ?? "—"))}</span>
      <span class="neighbor-score">${score.label} ${scoreText}</span>
      ${highBadge}
    </span>
  </button>`;
}

function cytoscapeStyles() {
  return [
    {
      selector: "node",
      style: {
        "background-color": "data(color)",
        "border-color": "data(borderColor)",
        "border-width": 2,
        "label": "data(label)",
        "font-size": 10,
        "font-family": "Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif",
        "font-weight": 700,
        "color": "#111111",
        "text-outline-color": "#ffffff",
        "text-outline-width": 2,
        "text-wrap": "none",
        "width": 22,
        "height": 22,
        "transition-property": "width height opacity border-width",
        "transition-duration": "240ms"
      }
    },
    {
      selector: "edge",
      style: {
        "curve-style": "bezier",
        "line-color": "data(edgeColor)",
        "width": 1.25,
        "opacity": .26,
        "transition-property": "opacity width line-color",
        "transition-duration": "220ms"
      }
    },
    {
      selector: "node:selected",
      style: {
        "border-width": 5,
        "width": 38,
        "height": 38,
        "font-size": 12,
        "font-weight": 850,
        "z-index": 9999
      }
    },
    {
      selector: ".focusNode",
      style: {
        "width": 44,
        "height": 44,
        "border-width": 6,
        "font-size": 12,
        "font-weight": 850,
        "opacity": 1,
        "z-index": 9999
      }
    },
    {
      selector: ".sharedFocusNode",
      style: {
        "background-color": "#ffffff",
        "pie-size": "100%",
        "pie-1-background-color": "data(primaryColor)",
        "pie-1-background-size": "50%",
        "pie-2-background-color": "data(secondaryColor)",
        "pie-2-background-size": "50%",
        "border-color": "#111827",
        "border-width": 5,
        "width": 46,
        "height": 46,
        "font-size": 12,
        "font-weight": 900,
        "color": "#000000",
        "opacity": 1,
        "z-index": 9999
      }
    },
    {
      selector: ".neighborNode",
      style: {
        "width": 29,
        "height": 29,
        "opacity": 1,
        "font-size": 11,
        "font-weight": 800,
        "color": "#000000",
        "z-index": 500
      }
    },
    {
      selector: ".sharedNeighborNode",
      style: {
        "background-color": "#ffffff",
        "pie-size": "100%",
        "pie-1-background-color": "data(primaryColor)",
        "pie-1-background-size": "50%",
        "pie-2-background-color": "data(secondaryColor)",
        "pie-2-background-size": "50%",
        "border-color": "#1f2937",
        "border-width": 2.5,
        "width": 32,
        "height": 32,
        "opacity": 1,
        "font-size": 11,
        "font-weight": 850,
        "color": "#000000",
        "z-index": 700
      }
    },
    {
      selector: ".dimNode",
      style: {
        "width": 8,
        "height": 8,
        "opacity": .34,
        "font-size": 0,
        "border-width": 1.2
      }
    },
    {
      selector: ".focusEdge",
      style: {
        "opacity": .9,
        "width": 2.8,
        "z-index": 500
      }
    },
    {
      selector: ".dimEdge",
      style: {
        "opacity": .035,
        "width": .5
      }
    },
    {
      selector: ".filtered",
      style: {
        "display": "none"
      }
    }
  ];
}

function destroyCy() {
  if (state.cy) {
    state.cy.destroy();
    state.cy = null;
  }
}

async function renderSingleSubsite(key, preserveSelection = false) {
  state.currentSubsite = key;
  state.activeSubsites = [key];
  state.crossMode = false;

  $("loading").classList.remove("hidden");
  const data = await loadSubsite(key);
  const cfg = APP_CONFIG.subsites[key];

  const nodeElements = [...data.connected].map(species => {
    const row = data.nodeMap.get(species) || { species };
    const cluster = row.cluster ?? "NA";
    return {
      data: {
        id: species,
        species,
        label: cleanSpeciesName(species),
        cluster: String(cluster),
        subsite: key,
        color: clusterColor(key, cluster),
        borderColor: cfg.subsiteColor,
        axis1: Number(row.Axis1),
        axis2: Number(row.Axis2)
      }
    };
  });

  const edgeElements = data.edges.map(e => ({
    data: {
      ...e,
      edgeColor: cfg.subsiteColor,
      subsite: key
    }
  }));

  destroyCy();

  state.cy = cytoscape({
    container: $("cy"),
    elements: [...nodeElements, ...edgeElements],
    style: cytoscapeStyles(),
    minZoom: .08,
    maxZoom: 4,
    wheelSensitivity: .22,
    selectionType: "single"
  });

  bindCyEvents();
  populateClusterFilter(key);
  populateSearch(key);
  renderClusterLegend([key]);
  renderSubsiteLegend([key]);
  updateHeader();
  applyLayout($("layoutMode").value);
  $("loading").classList.add("hidden");

  if (preserveSelection && state.selectedSpecies && data.connected.has(state.selectedSpecies)) {
    focusSpecies(key, state.selectedSpecies);
  } else {
    clearInfo();
  }
}

function bindCyEvents() {
  state.cy.on("tap", "node", evt => {
    const species = evt.target.data("species");
    const subsite = evt.target.data("subsite") || state.currentSubsite;

    if (state.crossMode) {
      // In cross mode the central focus node retains the selected species.
      if (evt.target.id() === "__focus__") return;
      // A shared neighbour is represented once with split cluster colours.
      // Use the primary subsite for its single-species detail panel rather than
      // trying to load a non-existent "shared" dataset.
      if (subsite === "shared") {
        showSpeciesInfo(state.activeSubsites[0], species);
      } else {
        showContextNodeInfo(subsite, species);
      }
    } else {
      focusSpecies(subsite, species);
    }
  });

  state.cy.on("tap", evt => {
    if (evt.target === state.cy) clearFocusClasses();
  });
}

function updateHeader() {
  const label = APP_CONFIG.subsites[state.currentSubsite].label;
  $("networkTitle").textContent = state.crossMode
    ? `${cleanSpeciesName(state.selectedSpecies)} — cross-subsite neighbourhood`
    : `${label} microbial network`;

  $("networkSubtitle").textContent = state.crossMode
    ? state.activeSubsites.map(k => APP_CONFIG.subsites[k].label).join(" + ")
    : "Click a node to reveal its direct ecological neighbourhood.";

  $("modeBadge").textContent = state.crossMode ? "Cross-subsite focus" : "Single subsite";
  $("subsiteLegendWrap").classList.toggle("hidden", !state.crossMode);

  document.querySelectorAll("#subsiteButtons button").forEach(btn => {
    btn.classList.toggle("active", btn.dataset.subsite === state.currentSubsite);
  });
}

function clearFocusClasses() {
  if (!state.cy || state.crossMode) return;
  state.cy.elements().removeClass("focusNode neighborNode dimNode focusEdge dimEdge");
}

function focusSpecies(subsite, species) {
  if (state.crossMode || subsite !== state.currentSubsite) {
    renderSingleSubsite(subsite).then(() => focusSpecies(subsite, species));
    return;
  }

  const cy = state.cy;
  const node = cy.getElementById(species);
  if (!node.length) return;

  state.selectedSpecies = species;
  clearFocusClasses();

  const neighbors = node.neighborhood("node");
  const directEdges = node.connectedEdges();

  node.addClass("focusNode");
  neighbors.addClass("neighborNode");
  directEdges.addClass("focusEdge");

  cy.nodes().difference(node.union(neighbors)).addClass("dimNode");
  cy.edges().difference(directEdges).addClass("dimEdge");

  arrangeFocus(node, neighbors);
  showSpeciesInfo(subsite, species);
}

function arrangeFocus(node, neighbors) {
  const cy = state.cy;
  const center = { x: cy.width() / 2, y: cy.height() / 2 };
  const n = Math.max(neighbors.length, 1);
  const radius = Math.min(245, Math.max(120, 48 + n * 7));

  node.animate({ position: center }, { duration: 350 });

  neighbors.forEach((nbr, i) => {
    const angle = (2 * Math.PI * i / n) - Math.PI / 2;
    nbr.animate({
      position: {
        x: center.x + radius * Math.cos(angle),
        y: center.y + radius * Math.sin(angle)
      }
    }, { duration: 380 });
  });

  const others = cy.nodes().difference(node.union(neighbors));
  const outerRadius = Math.max(cy.width(), cy.height()) * .48;
  const m = Math.max(others.length, 1);

  others.forEach((other, i) => {
    const angle = 2 * Math.PI * i / m;
    other.animate({
      position: {
        x: center.x + outerRadius * Math.cos(angle),
        y: center.y + outerRadius * Math.sin(angle)
      }
    }, { duration: 420 });
  });

  setTimeout(() => cy.fit(node.union(neighbors), 90), 430);
}

function applyLayout(mode) {
  if (!state.cy || state.crossMode) return;

  if (mode === "pca") {
    const scale = 24;
    state.cy.nodes().positions(node => ({
      x: state.cy.width() / 2 + (Number(node.data("axis1")) || 0) * scale,
      y: state.cy.height() / 2 + (Number(node.data("axis2")) || 0) * scale
    }));
    state.cy.fit(state.cy.elements(), 55);
  } else {
    state.cy.layout({
      name: "cose",
      animate: true,
      animationDuration: 650,
      fit: true,
      padding: 55,
      nodeRepulsion: 8500,
      idealEdgeLength: 80,
      edgeElasticity: 90,
      nestingFactor: 1.2,
      gravity: .22,
      numIter: 700,
      randomize: true
    }).run();
  }
}

function metadataFieldsForRow(subsite, row) {
  const cfg = APP_CONFIG.subsites[subsite];
  const preferred = ["cluster", ...cfg.scoreFields];
  const seen = new Set(preferred);

  const extras = Object.keys(row).filter(field => {
    if (seen.has(field) || field === "species" || field === "") return false;
    if (isHiddenMetadataField(field)) return false;
    const value = row[field];
    return value !== null && value !== undefined && value !== "";
  });

  return [...preferred, ...extras].filter(field => row[field] !== undefined);
}

async function showSpeciesInfo(subsite, species) {
  const data = await loadSubsite(subsite);
  const row = data.nodeMap.get(species);
  if (!row) return;

  $("infoEmpty").classList.add("hidden");
  $("clusterSummary").classList.add("hidden");
  $("infoContent").classList.remove("hidden");
  $("speciesTitle").textContent = cleanSpeciesName(species);

  const fields = metadataFieldsForRow(subsite, row);
  const focalScore = mainScoreInfo(subsite, row);
  const focalBadge = focalScore.high ? `<span class="high-score-badge focal-high-badge">${focalScore.highLabel}</span>` : "";
  $("speciesMeta").innerHTML = `
    <div class="meta-card">
      <div class="meta-heading"><span>${APP_CONFIG.subsites[subsite].label}</span>${focalBadge}</div>
      ${fields.map(field => `
        <div class="meta-row">
          <span>${titleCaseField(field)}</span>
          <span>${fmt(row[field])}</span>
        </div>
      `).join("")}
    </div>`;

  const neighbors = getNeighbors(subsite, species);
  $("connectionCount").textContent = `${neighbors.length} connection${neighbors.length === 1 ? "" : "s"}`;
  $("neighborList").innerHTML = neighbors
    .sort((a,b) => a.localeCompare(b))
    .map(n => neighborRowHtml(subsite, n))
    .join("");

  $("neighborList").querySelectorAll(".neighbor").forEach(btn => {
    btn.addEventListener("click", () => focusSpecies(subsite, btn.dataset.species));
  });

  await renderCrossButtons(species, subsite);
}

function showContextNodeInfo(subsite, species) {
  showSpeciesInfo(subsite, species);
}

function getNeighbors(subsite, species) {
  const edges = state.cache[subsite].edges;
  const set = new Set();
  edges.forEach(e => {
    if (e.source === species) set.add(e.target);
    if (e.target === species) set.add(e.source);
  });
  return [...set];
}

async function renderCrossButtons(species, primarySubsite) {
  const wrap = $("crossSubsiteButtons");
  wrap.innerHTML = "";

  const keys = Object.keys(APP_CONFIG.subsites).filter(k => k !== primarySubsite);
  let available = 0;

  keys.forEach(key => {
    const data = state.cache[key];
    if (!data || !data.nodeMap.has(species) || !data.connected.has(species)) return;
    available++;
    const btn = document.createElement("button");
    btn.className = "button";
    btn.textContent = `Add ${APP_CONFIG.subsites[key].label}`;
    btn.addEventListener("click", () => renderCrossSubsite(primarySubsite, key, species));
    wrap.appendChild(btn);
  });

  if (!available) {
    wrap.innerHTML = `<span class="muted small">No other connected subsite network contains this species.</span>`;
  }
}

async function renderCrossSubsite(primary, secondary, species, keepFilters = false) {
  state.currentSubsite = primary;
  state.activeSubsites = [primary, secondary];
  state.selectedSpecies = species;
  state.crossMode = true;

  const a = await loadSubsite(primary);
  const b = await loadSubsite(secondary);
  const aRow = a.nodeMap.get(species);
  const bRow = b.nodeMap.get(species);
  if (!aRow || !bRow) return;

  if (!keepFilters || !state.crossClusterFilters[primary]) state.crossClusterFilters[primary] = "all";
  if (!keepFilters || !state.crossClusterFilters[secondary]) state.crossClusterFilters[secondary] = "all";

  const aAllNeighbors = getNeighbors(primary, species);
  const bAllNeighbors = getNeighbors(secondary, species);

  const passesFilter = (subsite, nbr) => {
    const wanted = state.crossClusterFilters[subsite] || "all";
    if (wanted === "all") return true;
    return String(state.cache[subsite].nodeMap.get(nbr)?.cluster) === String(wanted);
  };

  const aVisible = new Set(aAllNeighbors.filter(n => passesFilter(primary, n)));
  const bVisible = new Set(bAllNeighbors.filter(n => passesFilter(secondary, n)));

  const sharedNeighbors = [...aVisible].filter(n => bVisible.has(n)).sort((x,y)=>x.localeCompare(y));
  const primaryOnly = [...aVisible].filter(n => !bVisible.has(n)).sort((x,y)=>x.localeCompare(y));
  const secondaryOnly = [...bVisible].filter(n => !aVisible.has(n)).sort((x,y)=>x.localeCompare(y));

  const centerX = $("cy").clientWidth / 2;
  const centerY = $("cy").clientHeight / 2 + 330;

  const elements = [{
    data: {
      id: "__focus__",
      species,
      label: cleanSpeciesName(species),
      subsite: "shared-focus",
      primarySubsite: primary,
      secondarySubsite: secondary,
      primaryCluster: String(aRow.cluster),
      secondaryCluster: String(bRow.cluster),
      primaryColor: clusterColor(primary, aRow.cluster),
      secondaryColor: clusterColor(secondary, bRow.cluster),
      color: "#ffffff",
      borderColor: "#111827"
    },
    position: { x: centerX, y: centerY },
    classes: "sharedFocusNode"
  }];

  const addExclusiveGroup = (subsite, neighbors, side) => {
    const cfg = APP_CONFIG.subsites[subsite];
    const data = state.cache[subsite];
    const n = neighbors.length;
    if (!n) return;

    const cols = Math.max(1, Math.ceil(Math.sqrt(n / 2)));
    const rows = Math.ceil(n / cols);
    const xGap = 125;
    const yGap = 66;
    const zoneCenterX = side === "left" ? centerX - 610 : centerX + 610;
    const zoneCenterY = centerY - 340;

    neighbors.forEach((nbr, i) => {
      const row = data.nodeMap.get(nbr) || {};
      const cluster = row.cluster ?? "NA";
      const col = i % cols;
      const r = Math.floor(i / cols);
      const id = `${subsite}::${nbr}`;

      elements.push({
        data: {
          id,
          species: nbr,
          label: cleanSpeciesName(nbr),
          subsite,
          cluster: String(cluster),
          color: clusterColor(subsite, cluster),
          borderColor: cfg.subsiteColor
        },
        position: {
          x: zoneCenterX + (col - (cols - 1) / 2) * xGap,
          y: zoneCenterY + (r - (rows - 1) / 2) * yGap
        },
        classes: "neighborNode"
      });

      elements.push({
        data: {
          id: `${subsite}::__focus__::${nbr}`,
          source: "__focus__",
          target: id,
          subsite,
          edgeColor: cfg.subsiteColor
        },
        classes: "focusEdge"
      });
    });
  };

  const addSharedGroup = neighbors => {
    const n = neighbors.length;
    if (!n) return;

    const cols = Math.max(1, Math.ceil(Math.sqrt(n * 1.6)));
    const rows = Math.ceil(n / cols);
    const xGap = 130;
    const yGap = 68;
    const zoneCenterX = centerX;
    const zoneCenterY = centerY - 350;

    neighbors.forEach((nbr, i) => {
      const pRow = state.cache[primary].nodeMap.get(nbr) || {};
      const sRow = state.cache[secondary].nodeMap.get(nbr) || {};
      const pCluster = pRow.cluster ?? "NA";
      const sCluster = sRow.cluster ?? "NA";
      const col = i % cols;
      const r = Math.floor(i / cols);
      const id = `shared::${nbr}`;

      elements.push({
        data: {
          id,
          species: nbr,
          label: cleanSpeciesName(nbr),
          subsite: "shared",
          primarySubsite: primary,
          secondarySubsite: secondary,
          primaryCluster: String(pCluster),
          secondaryCluster: String(sCluster),
          primaryColor: clusterColor(primary, pCluster),
          secondaryColor: clusterColor(secondary, sCluster),
          color: "#ffffff",
          borderColor: "#1f2937"
        },
        position: {
          x: zoneCenterX + (col - (cols - 1) / 2) * xGap,
          y: zoneCenterY + (r - (rows - 1) / 2) * yGap
        },
        classes: "sharedNeighborNode"
      });

      elements.push({
        data: {
          id: `${primary}::__focus__::shared::${nbr}`,
          source: "__focus__",
          target: id,
          subsite: primary,
          edgeColor: APP_CONFIG.subsites[primary].subsiteColor
        },
        classes: "focusEdge"
      });
      elements.push({
        data: {
          id: `${secondary}::__focus__::shared::${nbr}`,
          source: "__focus__",
          target: id,
          subsite: secondary,
          edgeColor: APP_CONFIG.subsites[secondary].subsiteColor
        },
        classes: "focusEdge"
      });
    });
  };

  addExclusiveGroup(primary, primaryOnly, "left");
  addSharedGroup(sharedNeighbors);
  addExclusiveGroup(secondary, secondaryOnly, "right");

  destroyCy();

  state.cy = cytoscape({
    container: $("cy"),
    elements,
    style: cytoscapeStyles(),
    layout: { name: "preset", fit: true, padding: 70 },
    minZoom: .12,
    maxZoom: 4,
    wheelSensitivity: .22
  });

  bindCyEvents();
  renderClusterLegend([primary, secondary]);
  renderSubsiteLegend([primary, secondary]);
  updateHeader();
  showCrossFocusInfo(
    primary,
    secondary,
    species,
    primaryOnly,
    sharedNeighbors,
    secondaryOnly,
    aAllNeighbors.length,
    bAllNeighbors.length
  );
}

function sharedNeighborRowHtml(primary, secondary, species) {
  const pRow = state.cache[primary]?.nodeMap.get(species) || {};
  const sRow = state.cache[secondary]?.nodeMap.get(species) || {};
  const pScore = mainScoreInfo(primary, pRow);
  const sScore = mainScoreInfo(secondary, sRow);
  const pText = Number.isFinite(pScore.value) ? pScore.value.toFixed(3) : "—";
  const sText = Number.isFinite(sScore.value) ? sScore.value.toFixed(3) : "—";
  return `<button class="neighbor" data-species="${escapeHtml(species)}" data-shared="true">
    <span class="neighbor-name">${escapeHtml(cleanSpeciesName(species))}</span>
    <span class="neighbor-details">
      <span class="neighbor-cluster">${APP_CONFIG.subsites[primary].label} C${escapeHtml(String(pRow.cluster ?? "—"))}</span>
      <span class="neighbor-score">${pScore.label} ${pText}</span>
      <span class="neighbor-cluster">${APP_CONFIG.subsites[secondary].label} C${escapeHtml(String(sRow.cluster ?? "—"))}</span>
      <span class="neighbor-score">${sScore.label} ${sText}</span>
    </span>
  </button>`;
}

function showCrossFocusInfo(primary, secondary, species, primaryOnly, sharedNeighbors, secondaryOnly, totalA, totalB) {
  $("infoEmpty").classList.add("hidden");
  $("clusterSummary").classList.add("hidden");
  $("infoContent").classList.remove("hidden");
  $("speciesTitle").textContent = cleanSpeciesName(species);

  const rows = [
    [primary, state.cache[primary].nodeMap.get(species)],
    [secondary, state.cache[secondary].nodeMap.get(species)]
  ];

  const metaCards = rows.map(([subsite, row]) => {
    const fields = metadataFieldsForRow(subsite, row);
    const score = mainScoreInfo(subsite, row);
    const badge = score.high ? `<span class="high-score-badge focal-high-badge">${score.highLabel}</span>` : "";
    return `
      <div class="meta-card">
        <div class="meta-heading"><span>${APP_CONFIG.subsites[subsite].label}</span>${badge}</div>
        ${fields.map(field => `
          <div class="meta-row"><span>${titleCaseField(field)}</span><span>${fmt(row[field])}</span></div>
        `).join("")}
      </div>`;
  }).join("");

  const clusterOptions = subsite => {
    const clusters = [...new Set(state.cache[subsite].nodeRows.map(r => r.cluster).filter(v => v !== null && v !== undefined))]
      .sort((a,b) => Number(a)-Number(b));
    return `<option value="all">All clusters</option>` + clusters.map(c =>
      `<option value="${escapeHtml(String(c))}" ${String(state.crossClusterFilters[subsite]) === String(c) ? "selected" : ""}>Cluster ${escapeHtml(String(c))}</option>`
    ).join("");
  };

  $("speciesMeta").innerHTML = `${metaCards}
    <div class="cross-cluster-controls">
      <div class="cross-filter-title">Filter direct neighbours by cluster</div>
      <div class="cross-filter-grid">
        <label><span>${APP_CONFIG.subsites[primary].label}</span><select id="crossClusterPrimary">${clusterOptions(primary)}</select></label>
        <label><span>${APP_CONFIG.subsites[secondary].label}</span><select id="crossClusterSecondary">${clusterOptions(secondary)}</select></label>
      </div>
    </div>`;

  $("crossClusterPrimary").addEventListener("change", e => {
    state.crossClusterFilters[primary] = e.target.value;
    renderCrossSubsite(primary, secondary, species, true);
  });
  $("crossClusterSecondary").addEventListener("change", e => {
    state.crossClusterFilters[secondary] = e.target.value;
    renderCrossSubsite(primary, secondary, species, true);
  });

  const visiblePrimary = primaryOnly.length + sharedNeighbors.length;
  const visibleSecondary = secondaryOnly.length + sharedNeighbors.length;
  $("connectionCount").textContent = `${visiblePrimary}/${totalA} ${APP_CONFIG.subsites[primary].label} + ${visibleSecondary}/${totalB} ${APP_CONFIG.subsites[secondary].label}`;

  const sharedSection = sharedNeighbors.length ? `
    <div class="neighbor-subsite-label">Common to both</div>
    ${sharedNeighbors.map(n => sharedNeighborRowHtml(primary, secondary, n)).join("")}` : "";
  const primarySection = primaryOnly.length ? `
    <div class="neighbor-subsite-label">${APP_CONFIG.subsites[primary].label} only</div>
    ${primaryOnly.map(n => neighborRowHtml(primary, n)).join("")}` : "";
  const secondarySection = secondaryOnly.length ? `
    <div class="neighbor-subsite-label">${APP_CONFIG.subsites[secondary].label} only</div>
    ${secondaryOnly.map(n => neighborRowHtml(secondary, n)).join("")}` : "";

  $("neighborList").innerHTML = sharedSection + primarySection + secondarySection;

  $("neighborList").querySelectorAll(".neighbor").forEach(btn => {
    btn.addEventListener("click", () => {
      const id = btn.dataset.species;
      if (btn.dataset.shared === "true") {
        // Shared species is shown once in the list; use the primary subsite for the single-species detail view.
        showSpeciesInfo(primary, id);
      } else if (primaryOnly.includes(id)) {
        showSpeciesInfo(primary, id);
      } else {
        showSpeciesInfo(secondary, id);
      }
    });
  });

  $("crossSubsiteButtons").innerHTML = "";
  Object.keys(APP_CONFIG.subsites)
    .filter(k => !state.activeSubsites.includes(k))
    .forEach(key => {
      const data = state.cache[key];
      if (!data.nodeMap.has(species) || !data.connected.has(species)) return;
      const btn = document.createElement("button");
      btn.className = "button";
      btn.textContent = `Switch second subsite to ${APP_CONFIG.subsites[key].label}`;
      btn.addEventListener("click", () => {
        state.crossClusterFilters[key] = "all";
        renderCrossSubsite(primary, key, species);
      });
      $("crossSubsiteButtons").appendChild(btn);
    });
}

function clearInfo() {
  state.selectedSpecies = null;
  $("infoContent").classList.add("hidden");
  const selectedCluster = $("clusterFilter")?.value || "all";
  if (!state.crossMode && selectedCluster !== "all") {
    renderClusterSummary(state.currentSubsite, selectedCluster);
  } else {
    $("clusterSummary").classList.add("hidden");
    $("infoEmpty").classList.remove("hidden");
  }
}

function renderClusterSummary(subsite, cluster) {
  const summary = $("clusterSummary");
  if (!summary) return;

  $("infoContent").classList.add("hidden");
  $("infoEmpty").classList.add("hidden");

  if (!subsite || cluster === "all") {
    summary.classList.add("hidden");
    $("infoEmpty").classList.remove("hidden");
    return;
  }

  const cfg = APP_CONFIG.subsites[subsite];
  const data = state.cache[subsite];
  const scoreField = subsite === "saliva" ? "sHACK" : "HAC";
  const highLabel = subsite === "saliva" ? "High sHACK" : "High HAC";

  const speciesRows = [...data.connected]
    .map(species => ({ species, row: data.nodeMap.get(species) }))
    .filter(x => x.row && String(x.row.cluster) === String(cluster))
    .map(x => ({
      species: x.species,
      score: Number(x.row[scoreField]),
      high: Number.isFinite(Number(x.row[scoreField])) && Number(x.row[scoreField]) >= 0.90
    }))
    .sort((a, b) => {
      const av = Number.isFinite(a.score) ? a.score : -Infinity;
      const bv = Number.isFinite(b.score) ? b.score : -Infinity;
      return bv - av || a.species.localeCompare(b.species);
    });

  const highCount = speciesRows.filter(x => x.high).length;
  const clusterColor = cfg.clusterColors[String(cluster)] || "#94a3b8";

  summary.innerHTML = `
    <div class="cluster-summary-head">
      <div>
        <div class="eyebrow">${escapeHtml(cfg.label)} cluster</div>
        <h2><span class="cluster-summary-dot" style="background:${clusterColor}"></span>Cluster ${escapeHtml(String(cluster))}</h2>
      </div>
    </div>
    <div class="cluster-summary-counts">
      <div><strong>${speciesRows.length}</strong><span>network species</span></div>
      <div><strong>${highCount}</strong><span>${highLabel} (≥0.90)</span></div>
    </div>
    <div class="cluster-summary-label">Species and ${scoreField}</div>
    <div class="cluster-species-list">
      ${speciesRows.length ? speciesRows.map(x => `
        <button class="cluster-species-row" data-species="${escapeHtml(x.species)}">
          <span class="cluster-species-name">${escapeHtml(cleanSpeciesName(x.species))}</span>
          <span class="cluster-species-score">${scoreField} ${Number.isFinite(x.score) ? x.score.toFixed(3) : "—"}</span>
          ${x.high ? `<span class="high-score-badge">${highLabel}</span>` : ""}
        </button>`).join("") : `<div class="muted small">No connected species are present in this cluster.</div>`}
    </div>`;

  summary.classList.remove("hidden");

  summary.querySelectorAll(".cluster-species-row").forEach(btn => {
    btn.addEventListener("click", () => focusSpecies(subsite, btn.dataset.species));
  });
}

function populateSearch(subsite) {
  const data = state.cache[subsite];
  const opts = [...data.connected]
    .sort((a,b) => a.localeCompare(b))
    .map(s => `<option value="${escapeHtml(s)}">${cleanSpeciesName(s)}</option>`)
    .join("");
  $("speciesOptions").innerHTML = opts;
  $("speciesSearch").value = "";
}

function populateClusterFilter(subsite) {
  const rows = state.cache[subsite].nodeRows;
  const clusters = [...new Set(rows.map(r => r.cluster).filter(v => v !== null && v !== undefined))]
    .sort((a,b) => Number(a)-Number(b));

  $("clusterFilter").innerHTML = `<option value="all">All clusters</option>` +
    clusters.map(c => `<option value="${c}">Cluster ${c}</option>`).join("");
}

function applyClusterFilter() {
  if (!state.cy || state.crossMode) return;
  const val = $("clusterFilter").value;
  state.cy.nodes().removeClass("filtered");
  state.cy.edges().removeClass("filtered");

  state.selectedSpecies = null;
  $("infoContent").classList.add("hidden");

  if (val === "all") {
    $("clusterSummary").classList.add("hidden");
    $("infoEmpty").classList.remove("hidden");
    state.cy.fit(state.cy.elements(), 55);
    return;
  }

  const hiddenNodes = state.cy.nodes().filter(n => String(n.data("cluster")) !== String(val));
  hiddenNodes.addClass("filtered");
  state.cy.edges().filter(e => hiddenNodes.contains(e.source()) || hiddenNodes.contains(e.target())).addClass("filtered");
  state.cy.fit(state.cy.nodes().not(".filtered"), 55);
  renderClusterSummary(state.currentSubsite, val);
}

function renderClusterLegend(subsites) {
  const wrap = $("clusterLegend");
  wrap.innerHTML = "";

  subsites.forEach(subsite => {
    const cfg = APP_CONFIG.subsites[subsite];
    Object.entries(cfg.clusterColors).forEach(([cluster, color]) => {
      const div = document.createElement("div");
      div.className = "legend-item";
      div.innerHTML = `<span class="legend-dot" style="background:${color}"></span>${cfg.label} C${cluster}`;
      wrap.appendChild(div);
    });
  });
}

function renderSubsiteLegend(subsites) {
  const wrap = $("subsiteLegend");
  wrap.innerHTML = "";
  subsites.forEach(subsite => {
    const cfg = APP_CONFIG.subsites[subsite];
    const div = document.createElement("div");
    div.className = "legend-item";
    div.innerHTML = `<span class="legend-dot" style="background:${cfg.subsiteColor}"></span>${cfg.label}`;
    wrap.appendChild(div);
  });
}

function renderSubsiteButtons() {
  const wrap = $("subsiteButtons");
  wrap.innerHTML = "";
  Object.entries(APP_CONFIG.subsites).forEach(([key, cfg]) => {
    const btn = document.createElement("button");
    btn.dataset.subsite = key;
    btn.textContent = cfg.label;
    btn.addEventListener("click", () => {
      state.selectedSpecies = null;
      $("layoutMode").value = "network";
      renderSingleSubsite(key);
    });
    wrap.appendChild(btn);
  });
}

function findSpecies() {
  const raw = $("speciesSearch").value.trim();
  if (!raw) return;

  const data = state.cache[state.currentSubsite];
  let species = raw;

  if (!data.connected.has(species)) {
    const norm = raw.toLowerCase().replaceAll(" ", "_");
    species = [...data.connected].find(s => s.toLowerCase() === norm) ||
              [...data.connected].find(s => s.toLowerCase().includes(norm));
  }

  if (species) focusSpecies(state.currentSubsite, species);
}

function numericMetricFields(subsite) {
  const cfg = APP_CONFIG.subsites[subsite];
  const rows = state.cache[subsite].nodeRows;

  const preferred = cfg.scoreFields.filter(f => rows.some(r => Number.isFinite(Number(r[f]))));

  // Automatically include future numeric fields (e.g. propensity) while excluding
  // PCA coordinates and redundant flag/module columns.
  const extras = [];
  Object.keys(rows[0] || {}).forEach(field => {
    if (preferred.includes(field) || field === "cluster" || isHiddenMetadataField(field)) return;
    const nums = rows.map(r => Number(r[field])).filter(Number.isFinite);
    if (nums.length >= Math.max(5, Math.floor(rows.length * .5))) extras.push(field);
  });

  return [...preferred, ...extras];
}

function openStats() {
  $("statsPanel").classList.remove("hidden");
  populateStatsSubsites();
  $("statsSubsite").value = state.currentSubsite;
  populateStatsMetrics(state.currentSubsite);
  renderStats();
  $("statsPanel").scrollIntoView({ behavior: "smooth", block: "start" });
}

function populateStatsSubsites() {
  $("statsSubsite").innerHTML = Object.entries(APP_CONFIG.subsites)
    .map(([k,cfg]) => `<option value="${k}">${cfg.label}</option>`)
    .join("");
}

function populateStatsMetrics(subsite) {
  const fields = numericMetricFields(subsite);
  $("statsMetric").innerHTML = fields
    .map(f => `<option value="${f}">${titleCaseField(f)}</option>`)
    .join("");
}

function median(values) {
  const a = [...values].sort((x,y) => x-y);
  if (!a.length) return NaN;
  const m = Math.floor(a.length/2);
  return a.length % 2 ? a[m] : (a[m-1]+a[m])/2;
}

function mean(values) {
  return values.reduce((a,b)=>a+b,0)/values.length;
}

function rankWithTies(values) {
  const indexed = values.map((v,i)=>({v,i})).sort((a,b)=>a.v-b.v);
  const ranks = Array(values.length);
  let i = 0;
  while (i < indexed.length) {
    let j = i + 1;
    while (j < indexed.length && indexed[j].v === indexed[i].v) j++;
    const avgRank = ((i + 1) + j) / 2;
    for (let k=i; k<j; k++) ranks[indexed[k].i] = avgRank;
    i = j;
  }
  return ranks;
}

function kruskalWallis(groups) {
  const flat = [];
  const groupIndex = [];
  groups.forEach((g, gi) => g.forEach(v => { flat.push(v); groupIndex.push(gi); }));
  const N = flat.length;
  if (N < 2 || groups.length < 2) return { H: NaN, df: 0, p: NaN };

  const ranks = rankWithTies(flat);
  const sums = Array(groups.length).fill(0);
  const ns = groups.map(g => g.length);
  ranks.forEach((r,i) => sums[groupIndex[i]] += r);

  let H = 12/(N*(N+1)) * sums.reduce((acc,R,i)=> acc + (R*R/ns[i]),0) - 3*(N+1);

  // Tie correction
  const counts = new Map();
  flat.forEach(v => counts.set(v, (counts.get(v)||0)+1));
  const tieTerm = [...counts.values()].reduce((a,t)=> a + (t*t*t - t),0);
  const correction = 1 - tieTerm/(N*N*N - N);
  if (correction > 0) H /= correction;

  const df = groups.length - 1;
  const p = 1 - jStat.chisquare.cdf(H, df);
  return { H, df, p };
}

function renderStats() {
  const subsite = $("statsSubsite").value;
  const metric = $("statsMetric").value;
  if (!subsite || !metric) return;

  const cfg = APP_CONFIG.subsites[subsite];
  const rows = state.cache[subsite].nodeRows
    .filter(r => r.cluster !== null && r.cluster !== undefined && Number.isFinite(Number(r[metric])));

  const clusters = [...new Set(rows.map(r => String(r.cluster)))]
    .sort((a,b)=>Number(a)-Number(b));

  const traces = clusters.map(cluster => {
    const subset = rows.filter(r => String(r.cluster) === cluster);
    return {
      type: "box",
      name: `Cluster ${cluster}`,
      y: subset.map(r => Number(r[metric])),
      boxpoints: "all",
      jitter: .35,
      pointpos: 0,
      marker: { color: clusterColor(subsite, cluster), size: 5, opacity: .72 },
      line: { color: clusterColor(subsite, cluster), width: 1.5 },
      fillcolor: hexToRgba(clusterColor(subsite, cluster), .18),
      customdata: subset.map(r => speciesFromRow(r)),
      hovertemplate: "%{customdata}<br>" + metric + ": %{y:.3f}<extra></extra>"
    };
  });

  Plotly.react("statsPlot", traces, {
    margin: { l: 55, r: 20, t: 30, b: 55 },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(0,0,0,0)",
    showlegend: false,
    xaxis: { title: "Cluster", showgrid: false, zeroline: false },
    yaxis: { title: titleCaseField(metric), gridcolor: "#e8edf4", zeroline: false },
    font: { family: "Inter, system-ui, sans-serif", color: "#475569" },
    hoverlabel: { bgcolor: "#111827", font: { color: "#fff" } }
  }, { responsive: true, displaylogo: false });

  const plot = $("statsPlot");
  if (plot.removeAllListeners) plot.removeAllListeners("plotly_click");
  plot.on("plotly_click", ev => {
    const species = ev?.points?.[0]?.customdata;
    if (!species) return;
    renderSingleSubsite(subsite).then(() => focusSpecies(subsite, species));
    window.scrollTo({ top: 0, behavior: "smooth" });
  });

  $("statsTitle").textContent = `${titleCaseField(metric)} across ${cfg.label} clusters`;

  const groups = clusters.map(c => rows.filter(r => String(r.cluster) === c).map(r => Number(r[metric])));
  const kw = kruskalWallis(groups);
  $("kwResult").textContent = Number.isFinite(kw.p)
    ? `Kruskal–Wallis H = ${kw.H.toFixed(2)}, df = ${kw.df}, p ${kw.p < .001 ? "< 0.001" : "= " + kw.p.toFixed(4)}`
    : "Kruskal–Wallis test unavailable";

  $("summaryTable").querySelector("tbody").innerHTML = clusters.map((c,i) => {
    const g = groups[i];
    return `<tr>
      <td>Cluster ${c}</td>
      <td>${g.length}</td>
      <td>${mean(g).toFixed(3)}</td>
      <td>${median(g).toFixed(3)}</td>
    </tr>`;
  }).join("");
}

function hexToRgba(hex, a) {
  const h = hex.replace("#","");
  const bigint = parseInt(h,16);
  const r = (bigint>>16)&255, g=(bigint>>8)&255, b=bigint&255;
  return `rgba(${r},${g},${b},${a})`;
}

function escapeHtml(str) {
  return String(str)
    .replaceAll("&","&amp;")
    .replaceAll("<","&lt;")
    .replaceAll(">","&gt;")
    .replaceAll('"',"&quot;")
    .replaceAll("'","&#039;");
}

async function init() {
  renderSubsiteButtons();
  $("loading").classList.remove("hidden");

  try {
    await preloadAll();

    $("clusterFilter").addEventListener("change", applyClusterFilter);
    $("layoutMode").addEventListener("change", e => applyLayout(e.target.value));
    $("resetButton").addEventListener("click", () => {
      state.selectedSpecies = null;
      $("layoutMode").value = "network";
      renderSingleSubsite(state.currentSubsite);
    });
    $("searchButton").addEventListener("click", findSpecies);
    $("speciesSearch").addEventListener("keydown", e => { if (e.key === "Enter") findSpecies(); });
    $("closeFocus").addEventListener("click", () => {
      if (state.crossMode) renderSingleSubsite(state.currentSubsite);
      else {
        clearInfo();
        clearFocusClasses();
        applyLayout($("layoutMode").value);
      }
    });

    $("statsToggle").addEventListener("click", openStats);
    $("statsClose").addEventListener("click", () => $("statsPanel").classList.add("hidden"));
    $("statsSubsite").addEventListener("change", e => {
      populateStatsMetrics(e.target.value);
      renderStats();
    });
    $("statsMetric").addEventListener("change", renderStats);

    await renderSingleSubsite(APP_CONFIG.defaultSubsite);
  } catch (err) {
    console.error(err);
    $("loading").textContent = "Could not load the CSV files. Open this project through a web server (for example GitHub Pages), not by double-clicking index.html.";
    $("loading").classList.remove("hidden");
  }
}

document.addEventListener("DOMContentLoaded", init);
