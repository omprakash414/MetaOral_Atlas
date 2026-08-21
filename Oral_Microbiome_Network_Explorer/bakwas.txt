# Oral Microbiome Network Explorer

A browser-based interactive explorer for subsite-specific oral microbiome networks.

## What it does

- Switch among Saliva, Supragingival, Subgingival and Tongue–tonsil networks.
- Colour nodes by the **cluster/module number within each subsite**.
- Search for a species and focus the graph on its direct neighbours.
- Selecting a specific cluster shows a right-panel cluster summary with species-level sHACK/HAC scores and the number of High sHACK/High HAC taxa (score ≥ 0.90).
- On node selection, direct neighbours remain prominent while unrelated nodes shrink, fade and move outward.
- Show the selected species' available scores and metadata.
- Detect species shared across subsites.
- Add a second subsite around a shared focal species and display the direct neighbourhoods from both subsites separately.
- Compare score distributions across clusters using interactive boxplots and a Kruskal–Wallis test.
- Click a point in a statistics plot to return to that species in the network.
- The network view uses the **Compound Spring Embedder (CoSE)** layout for the edge-based network display.
- The **Direct connections** list supports horizontal scrolling so long species names remain visible.
- In two-subsite comparison mode, direct neighbours are separated into **primary-only**, **common to both**, and **secondary-only** groups.
- A species directly connected in both active subsites is shown **once** with a 50/50 split node colour corresponding to its cluster in each subsite.
- Network species labels are kept on a **single line**.
- Switch between a force-directed network layout and the Axis1/Axis2 association-space layout.

## Repository structure

```text
Oral_Microbiome_Network_Explorer/
├── index.html
├── .nojekyll
├── README.md
├── css/
│   └── style.css
├── js/
│   ├── config.js
│   └── app.js
└── data/
    ├── saliva/
    │   ├── S8_1Saliva_PCAdf.csv
    │   └── S7_1saliva_ControlCohort_RemNetwork_melt_filt_0.0001.csv
    ├── supragingival/
    │   ├── S8_1Supragingival_PCAdf.csv
    │   └── S7_4supragingival_ControlCohort_RemNetwork_melt_filt.csv
    ├── subgingival/
    │   ├── S8_1Subgingival_PCAdf(1).csv
    │   └── S7_5subgingival_ControlCohort_RemNetwork_melt_filt.csv
    └── tongue_tonsil/
        ├── S8_1Tongue_PCAdf(1).csv
        └── S7_6tongue_tonsil_ControlCohort_RemNetwork_melt_filt.csv
```

## Important data conventions

### Node/value files
The first column (`Unnamed: 0` in the current files) is treated as the species identifier.

The application deliberately does **not** display the redundant `High_HAC`, `High_sHACK` or `Module*` flags in the species panel. It shows the actual score(s) and `cluster` number instead.

Current default score fields are:

- **Saliva:** `CS`, `HS`, `SS`, `sHACK`
- **Supragingival:** `HAC`
- **Subgingival:** `HAC`
- **Tongue–tonsil:** `HAC`

If additional numeric node properties such as propensity are later added to these CSVs, the application will automatically make them available to the metadata/statistics layer unless their field names match one of the explicitly hidden technical columns.

### Edge files
The current edge files use:

- `Node1`
- `Node2`
- `Weight`

The current `Weight` values can remain as they are.

## Run locally

Do not open `index.html` by double-clicking it because browsers can block local CSV loading.

From the repository directory, run for example:

```bash
python3 -m http.server 8000
```

Then open:

```text
http://localhost:8000
```

## Deploy with GitHub Pages

1. Create a GitHub repository, for example `Oral_Microbiome_Network_Explorer`.
2. Upload/push this complete folder structure.
3. On GitHub, open **Settings → Pages**.
4. Under **Build and deployment**, select **Deploy from a branch**.
5. Choose the `main` branch and `/ (root)`.
6. Save.
7. GitHub will publish a URL in the form:

```text
https://YOUR_GITHUB_USERNAME.github.io/Oral_Microbiome_Network_Explorer/
```

## Add your own domain later

A custom domain can be connected to GitHub Pages after the explorer is working. Configure the domain only after the GitHub Pages URL works correctly.

## Libraries loaded in the browser

The application uses CDN-hosted versions of:

- Cytoscape.js
- Papa Parse
- Plotly.js
- jStat

No Python, R, database, Node.js server or backend is required for visitors to use the explorer.


### Species identifier handling

The supplied `*_PCAdf.csv` files store species names in the unnamed first CSV column. The application reads this column directly and joins it to `Node1`/`Node2` in the edge files. This is required for correct cluster colours and the right-hand species information panel.

## Interaction refinements

- `sHACK >= 0.90` is labelled **High sHACK** for Saliva; `HAC >= 0.90` is labelled **High HAC** for the other subsites.
- Direct-connection rows show the neighbour species, cluster, sHACK/HAC score and high-score badge inside the existing scroll area.
- Cross-subsite neighbourhood mode has independent cluster filters for each of the two active subsites.
- Non-neighbour nodes remain visible as small, faint nodes during focus mode.
- Species labels use darker, more readable text while retaining the original site theme and layout.
