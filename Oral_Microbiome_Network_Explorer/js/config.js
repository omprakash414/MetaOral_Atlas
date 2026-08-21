const APP_CONFIG = {
  defaultSubsite: "saliva",

  subsites: {
    saliva: {
      label: "Saliva",
      nodeFile: "data/saliva/S8_1Saliva_PCAdf.csv",
      edgeFile: "data/saliva/S7_1saliva_ControlCohort_RemNetwork_melt_filt_0.0001.csv",
      scoreFields: ["CS", "HS", "SS", "sHACK"],
      clusterColors: {
        "1": "#1f77b4",
        "2": "#ff7f0e",
        "3": "#2ca02c",
        "4": "#d62728",
        "5": "#9467bd"
      },
      subsiteColor: "#1d4ed8"
    },

    supragingival: {
      label: "Supragingival",
      nodeFile: "data/supragingival/S8_1Supragingival_PCAdf.csv",
      edgeFile: "data/supragingival/S7_4supragingival_ControlCohort_RemNetwork_melt_filt.csv",
      scoreFields: ["HAC"],
      clusterColors: {
        "1": "#17becf",
        "2": "#bcbd22",
        "3": "#8c564b"
      },
      subsiteColor: "#c2410c"
    },

    subgingival: {
      label: "Subgingival",
      nodeFile: "data/subgingival/S8_1Subgingival_PCAdf(1).csv",
      edgeFile: "data/subgingival/S7_5subgingival_ControlCohort_RemNetwork_melt_filt.csv",
      scoreFields: ["HAC"],
      clusterColors: {
        "1": "#e377c2",
        "2": "#7f7f7f"
      },
      subsiteColor: "#15803d"
    },

    tongue_tonsil: {
      label: "Tongue–tonsil",
      nodeFile: "data/tongue_tonsil/S8_1Tongue_PCAdf(1).csv",
      edgeFile: "data/tongue_tonsil/S7_6tongue_tonsil_ControlCohort_RemNetwork_melt_filt.csv",
      scoreFields: ["HAC"],
      clusterColors: {
        "1": "#003f5c",
        "2": "#58508d",
        "3": "#bc5090",
        "4": "#ff6361",
        "5": "#ffa600",
        "6": "#2f4b7c"
      },
      subsiteColor: "#be185d"
    }
  },

  // These fields are intentionally hidden from the species information panel.
  // High_HAC / High_sHACK and Module* are redundant with the actual score
  // and cluster number, as requested.
  hiddenMetadataPatterns: [
    /^$/,
    /^Unnamed:/i,
    /^Axis[123]$/i,
    /^c[xy]$/i,
    /^cluster_factor$/i,
    /^High_/i,
    /^Module\d+$/i
  ]
};
