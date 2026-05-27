const RESTRICTION_SITES = {
  HindIII: "AAGCTT",
  BamHI: "GGATCC",
  XhoI: "CTCGAG",
  SalI: "GTCGAC",
  KpnI: "GGTACC",
  BglII: "AGATCT",
};

const CODON_TABLE = {
  TTT: "F", TTC: "F", TTA: "L", TTG: "L", CTT: "L", CTC: "L", CTA: "L", CTG: "L",
  ATT: "I", ATC: "I", ATA: "I", ATG: "M",
  GTT: "V", GTC: "V", GTA: "V", GTG: "V",
  TCT: "S", TCC: "S", TCA: "S", TCG: "S", AGT: "S", AGC: "S",
  CCT: "P", CCC: "P", CCA: "P", CCG: "P",
  ACT: "T", ACC: "T", ACA: "T", ACG: "T",
  GCT: "A", GCC: "A", GCA: "A", GCG: "A",
  TAT: "Y", TAC: "Y",
  TAA: "*", TAG: "*", TGA: "*",
  CAT: "H", CAC: "H", CAA: "Q", CAG: "Q",
  AAT: "N", AAC: "N", AAA: "K", AAG: "K",
  GAT: "D", GAC: "D", GAA: "E", GAG: "E",
  TGT: "C", TGC: "C", TGG: "W",
  CGT: "R", CGC: "R", CGA: "R", CGG: "R", AGA: "R", AGG: "R",
  GGT: "G", GGC: "G", GGA: "G", GGG: "G",
};

const OPTIMAL_CODONS = {
  A: "GCC", C: "TGC", D: "GAC", E: "GAG", F: "TTC",
  G: "GGC", H: "CAC", I: "ATC", K: "AAG", L: "CTG",
  M: "ATG", N: "AAC", P: "CCC", Q: "CAG", R: "AGG",
  S: "AGC", T: "ACC", V: "GTG", W: "TGG", Y: "TAC",
};

const HOST_CONFIG = {
  HEK: {
    label: "HEK",
    vectors: ["p3.4", "P3.3", "P3.4 C-term m2aFc"],
    signalPeptides: [
      { id: "h7", label: "H7", dna: "ATGGAATTCGGCCTGTCCTGGGTGTTCCTGGTGGCCCTGTTTAGAGGCGTGCAGTGC", aa: "MEFGLSWVFLVALFRGVQC" },
      { id: "gl", label: "GL", dna: "ATGGGAGTGAAGGTGCTGTTCGCCCTGATCTGTATCGCCGTGGCCGAGGCC", aa: "MGVKVLFALICIAVAEA" },
      { id: "hsa", label: "HSA", dna: "ATGAAGTGGGTCACCTTTATCTCCCTGCTGTTCCTGTTCTCCTCCGCCTACTCT", aa: "MKWVTFISLLFLFSSAYS" },
      { id: "igh", label: "IgH", dna: "ATGGGCTGGTCCTGCATCATCCTGTTTCTGGTGGCCACAGCCACAGGCGTGCACTCT", aa: "MGWSCIILFLVATATGVHS" },
      { id: "secrecon", label: "Secrecon", dna: "ATGTGGTGGAGGCTGTGGTGGCTGTTACTCCTGCTGCTGCTGCTGTGGCCCATGGTGTGGGCA", aa: "MWWRLWWLLLLLLLLWPMVWA" },
    ],
  },
  CHO: {
    label: "CHO",
    vectors: ["p3.4", "P3.3", "P3.4 C-term m2aFc"],
    signalPeptides: [
      { id: "h7", label: "H7", dna: "ATGGAATTCGGCCTGTCCTGGGTGTTCCTGGTGGCCCTGTTTAGAGGCGTGCAGTGC", aa: "MEFGLSWVFLVALFRGVQC" },
      { id: "gl", label: "GL", dna: "ATGGGAGTGAAGGTGCTGTTCGCCCTGATCTGTATCGCCGTGGCCGAGGCC", aa: "MGVKVLFALICIAVAEA" },
      { id: "hsa", label: "HSA", dna: "ATGAAGTGGGTCACCTTTATCTCCCTGCTGTTCCTGTTCTCCTCCGCCTACTCT", aa: "MKWVTFISLLFLFSSAYS" },
      { id: "igh", label: "IgH", dna: "ATGGGCTGGTCCTGCATCATCCTGTTTCTGGTGGCCACAGCCACAGGCGTGCACTCT", aa: "MGWSCIILFLVATATGVHS" },
      { id: "secrecon", label: "Secrecon", dna: "ATGTGGTGGAGGCTGTGGTGGCTGTTACTCCTGCTGCTGCTGCTGTGGCCCATGGTGTGGGCA", aa: "MWWRLWWLLLLLLLLWPMVWA" },
    ],
  },
  Baculovirus: {
    label: "Baculovirus (insect cell)",
    vectors: [
      "BacPAK9", "Bachis9", "BacGH9", "BacHG9", "Bachis9(HM)", "Bachis9(SP12)",
      "AceBac9", "AceBachis9", "AceBacGH9", "AceBacHG9", "AceBachis9(HM)", "AceBachis9(SP12)",
    ],
    signalPeptides: [
      { id: "hm", label: "H-M", dna: "ATGGTCTCTGCGATCGTCCTCTACGTGCTGCTGGCTGCCGCAGCTCACAGCGCTTTCGCT", aa: "MVSAIVLYVLLAAAAHSAFA" },
      { id: "sp12", label: "SP12", dna: "ATGCGCGTGCTGGTGCTGCTGGCGTGCCTGGCGGCGGCGAGCAACGCG", aa: "MRVLVLLACLAAASNA" },
    ],
  },
};

const KOZAK_SEQ = "GCCACC";
const HIS6_DNA = "GGATCTCATCATCATCACCATCAC";
const HIS6_AA = "GSHHHHHH";
const STOP_CODON = "TAA";
const EXTRA_5 = "GCCA";
const EXTRA_3 = "GCCA";
const DNA_COMPLEMENT = { A: "T", T: "A", C: "G", G: "C", N: "N" };

let proteinState = null;
let finalState = null;
let activeTab = "report";

const el = (id) => document.getElementById(id);

function cleanProtein(raw) {
  return raw
    .split(/\r?\n/)
    .filter((line) => !line.trim().startsWith(">"))
    .join("")
    .replace(/[^A-Za-z*]/g, "")
    .replace(/\*/g, "")
    .toUpperCase();
}

function cleanDna(raw) {
  return raw
    .split(/\r?\n/)
    .filter((line) => !line.trim().startsWith(">"))
    .join("")
    .replace(/[^ATCGatcgNn]/g, "")
    .toUpperCase();
}

function translateDna(dna) {
  const aa = [];
  for (let i = 0; i <= dna.length - 3; i += 3) {
    const residue = CODON_TABLE[dna.slice(i, i + 3)] || "?";
    if (residue === "*") break;
    aa.push(residue);
  }
  return aa.join("");
}

function reverseTranslate(protein) {
  return [...protein].map((aa) => OPTIMAL_CODONS[aa] || "NNN").join("");
}

function stripTerminalStop(dna) {
  const stop = dna.slice(-3);
  return ["TAA", "TAG", "TGA"].includes(stop) ? dna.slice(0, -3) : dna;
}

function reverseTranslateDesignedProtein(proteinState) {
  if (proteinState.config.host !== "Baculovirus" && proteinState.config.signalStrategy === "heterologous") {
    const sp = proteinState.config.signalPeptide;
    if (proteinState.designedProtein.startsWith(sp.aa)) {
      return `${sp.dna}${reverseTranslate(proteinState.designedProtein.slice(sp.aa.length))}`;
    }
  }
  return reverseTranslate(proteinState.targetProtein);
}

function formatSeq(seq, width = 105) {
  const lines = [];
  for (let i = 0; i < seq.length; i += width) {
    lines.push(seq.slice(i, i + width));
  }
  return lines.join("\n");
}

function fasta(name, seq, width = 105) {
  return `>${name}\n${formatSeq(seq, width)}`;
}

function reverseComplement(seq) {
  return seq
    .split("")
    .reverse()
    .map((base) => DNA_COMPLEMENT[base] || "N")
    .join("");
}

function findRestrictionSites(seq) {
  return Object.entries(RESTRICTION_SITES).map(([enzyme, site]) => {
    const positions = [];
    let start = 0;
    while (true) {
      const found = seq.indexOf(site, start);
      if (found === -1) break;
      positions.push(found + 1);
      start = found + 1;
    }
    return { enzyme, site, positions };
  });
}

function gcPercent(seq) {
  if (!seq.length) return 0;
  return ([...seq].filter((base) => base === "G" || base === "C").length / seq.length) * 100;
}

function wallaceTm(seq) {
  const counts = [...seq].reduce((acc, base) => {
    acc[base] = (acc[base] || 0) + 1;
    return acc;
  }, {});
  return 2 * ((counts.A || 0) + (counts.T || 0)) + 4 * ((counts.G || 0) + (counts.C || 0));
}

function chooseAnnealLength(seq, fromEnd = false) {
  const minLen = Math.min(18, seq.length);
  const maxLen = Math.min(24, seq.length);
  let best = minLen;
  let bestScore = Number.POSITIVE_INFINITY;
  for (let len = minLen; len <= maxLen; len += 1) {
    const segment = fromEnd ? seq.slice(-len) : seq.slice(0, len);
    const score = Math.abs(wallaceTm(segment) - 60) + Math.abs(gcPercent(segment) - 50) / 3;
    if (score < bestScore) {
      best = len;
      bestScore = score;
    }
  }
  return best;
}

function currentConfig() {
  const host = el("hostSelect").value;
  const vector = el("vectorSelect").value;
  const signalStrategy = el("signalStrategySelect").value;
  const signalPeptideId = el("signalPeptideSelect").value;
  const signalPeptide = HOST_CONFIG[host].signalPeptides.find((item) => item.id === signalPeptideId) || HOST_CONFIG[host].signalPeptides[0];
  return {
    host,
    vector,
    signalStrategy,
    signalPeptide,
    nativeSignalAa: cleanProtein(el("nativeSignalAa").value),
    cleavageProbability: el("cleavageProbability").value.trim(),
    nativeSignalLabel: el("nativeSignalLabel").value.trim(),
    tmStart: Number(el("tmStart").value),
    tmEnd: Number(el("tmEnd").value),
    infusionEnabled: el("infusionToggle").checked,
    forwardOverlap: cleanDna(el("forwardOverlap").value),
    reverseOverlap: cleanDna(el("reverseOverlap").value),
    codonTool: el("codonToolSelect").value,
    hostOrganism: el("hostOrganismSelect").value,
    accession: el("accessionInput").value.trim(),
    purpose: el("purposeInput").value.trim(),
    notes: el("notesInput").value.trim(),
  };
}

function vectorEncodedSignalPeptide(config) {
  if (config.host !== "Baculovirus") return null;
  if (config.vector.includes("(HM)")) return HOST_CONFIG.Baculovirus.signalPeptides.find((sp) => sp.id === "hm");
  if (config.vector.includes("(SP12)")) return HOST_CONFIG.Baculovirus.signalPeptides.find((sp) => sp.id === "sp12");
  return null;
}

function designProtein(inputProtein, config) {
  let designed = inputProtein;
  let targetProtein = inputProtein;
  const notes = [];
  let nativeMatched = false;

  if (config.nativeSignalAa && designed.startsWith(config.nativeSignalAa)) {
    nativeMatched = true;
    notes.push(`Native signal peptide: ${config.nativeSignalLabel || config.nativeSignalAa}`);
  }

  if (config.signalStrategy === "none") {
    if (nativeMatched) {
      const mature = designed.slice(config.nativeSignalAa.length);
      designed = `M${mature}`;
      targetProtein = designed;
      notes.push("Signal peptide removed for cytosol expression");
    } else {
      notes.push("No signal peptide added: cytosol expression");
    }
  }

  if (config.signalStrategy === "native") {
    targetProtein = designed;
    notes.push(nativeMatched ? "Native signal peptide retained" : "Native signal peptide strategy selected");
  }

  if (config.signalStrategy === "heterologous") {
    const mature = nativeMatched ? designed.slice(config.nativeSignalAa.length) : designed.replace(/^M/, "");
    designed = `${config.signalPeptide.aa}${mature}`;
    targetProtein = designed;
    notes.push(`${config.signalPeptide.label} signal peptide inserted`);
  }

  const vectorSp = vectorEncodedSignalPeptide(config);
  if (config.host === "Baculovirus") {
    const mature = nativeMatched ? inputProtein.slice(config.nativeSignalAa.length) : inputProtein.replace(/^M/, "");
    targetProtein = `M${mature}`;
    if (vectorSp) {
      designed = `${vectorSp.aa}GS${mature}`;
      notes.push(`${config.vector} vector-encoded ${vectorSp.label} signal peptide and GS linker noted; SP is excluded from insert`);
    } else {
      designed = targetProtein;
      notes.push("Baculovirus insert excludes signal peptide and vector C-terminal His tag");
    }
  }

  if (config.tmStart || config.tmEnd) {
    if (!config.tmStart || !config.tmEnd || config.tmStart < 1 || config.tmEnd < config.tmStart || config.tmEnd > designed.length) {
      throw new Error(`TM domain 범위는 1-${designed.length} 안에서 시작 <= 끝이어야 합니다.`);
    }
    const removed = designed.slice(config.tmStart - 1, config.tmEnd);
    designed = `${designed.slice(0, config.tmStart - 1)}${designed.slice(config.tmEnd)}`;
    targetProtein = `${targetProtein.slice(0, Math.min(config.tmStart - 1, targetProtein.length))}${targetProtein.slice(Math.min(config.tmEnd, targetProtein.length))}`;
    notes.push(`TM domain removed: AA ${config.tmStart}-${config.tmEnd} (${removed})`);
  }

  return { inputProtein, designedProtein: designed, targetProtein, notes, config };
}

function buildConstruct(cdsNoStop, config) {
  if (config.host === "Baculovirus") {
    return `${RESTRICTION_SITES.BamHI}${cdsNoStop}${RESTRICTION_SITES.XhoI}`;
  }
  return `${RESTRICTION_SITES.HindIII}${KOZAK_SEQ}${cdsNoStop}${HIS6_DNA}${STOP_CODON}${RESTRICTION_SITES.BamHI}`;
}

function annotate(cdsNoStop, config) {
  const parts = config.host === "Baculovirus"
    ? [
      ["BamHI", RESTRICTION_SITES.BamHI],
      ["Target DNA", cdsNoStop],
      ["XhoI", RESTRICTION_SITES.XhoI],
    ]
    : [
      ["HindIII", RESTRICTION_SITES.HindIII],
      ["Kozak", KOZAK_SEQ],
      ["SP + Target DNA", cdsNoStop],
      [`6xHis tag (${HIS6_AA})`, HIS6_DNA],
      ["Stop codon", STOP_CODON],
      ["BamHI", RESTRICTION_SITES.BamHI],
    ];
  let pos = 1;
  return parts.map(([label, seq]) => {
    const row = { label, start: pos, end: pos + seq.length - 1, seq };
    pos += seq.length;
    return row;
  });
}

function makePrimers(cdsNoStop, config) {
  const forwardAnneal = cdsNoStop.slice(0, chooseAnnealLength(cdsNoStop, false));
  const reverseAnnealCoding = cdsNoStop.slice(-chooseAnnealLength(cdsNoStop, true));
  const forwardTail = config.infusionEnabled
    ? `${config.forwardOverlap}${config.host === "Baculovirus" ? "" : KOZAK_SEQ}`
    : config.host === "Baculovirus"
      ? RESTRICTION_SITES.BamHI
      : `${RESTRICTION_SITES.HindIII}${KOZAK_SEQ}`;
  const reverseTemplate = config.infusionEnabled
    ? config.host === "Baculovirus"
      ? reverseAnnealCoding
      : `${reverseAnnealCoding}${HIS6_DNA}${STOP_CODON}`
    : config.host === "Baculovirus"
      ? `${reverseAnnealCoding}${RESTRICTION_SITES.XhoI}`
      : `${reverseAnnealCoding}${HIS6_DNA}${STOP_CODON}${RESTRICTION_SITES.BamHI}`;
  const reverseTail = config.infusionEnabled ? config.reverseOverlap : "";
  const reverseAnneal = reverseComplement(reverseAnnealCoding);

  return {
    mode: config.infusionEnabled ? "Infusion" : "Restriction enzyme",
    forward: {
      name: `${config.host}_${config.vector}_${config.infusionEnabled ? "Infusion" : config.host === "Baculovirus" ? "BamHI" : "HindIII"}_F`,
      sequence: `${forwardTail}${forwardAnneal}`,
      tail: forwardTail,
      anneal: forwardAnneal,
      tm: wallaceTm(forwardAnneal),
      gc: gcPercent(forwardAnneal),
    },
    reverse: {
      name: `${config.host}_${config.vector}_${config.infusionEnabled ? "Infusion" : config.host === "Baculovirus" ? "XhoI" : "BamHI"}_R`,
      sequence: `${reverseTail}${reverseComplement(reverseTemplate)}`,
      tail: reverseTail || reverseComplement(reverseTemplate).slice(0, reverseComplement(reverseTemplate).length - reverseAnneal.length),
      anneal: reverseAnneal,
      tm: wallaceTm(reverseAnneal),
      gc: gcPercent(reverseAnneal),
    },
  };
}

function finalAnalyze() {
  if (!proteinState) throw new Error("먼저 아미노산 서열을 적용하세요.");
  const config = currentConfig();
  const dnaInput = cleanDna(el("optimizedDnaInput").value);
  const cdsNoStop = stripTerminalStop(dnaInput || proteinState.reverseTranslatedDna || "");
  if (!cdsNoStop) throw new Error("Reverse translation 또는 codon optimization DNA가 필요합니다.");
  if (!cdsNoStop.startsWith("ATG")) throw new Error("최종 CDS는 ATG로 시작해야 합니다.");
  if (cdsNoStop.length % 3 !== 0) throw new Error("최종 CDS 길이는 3의 배수여야 합니다.");

  const translated = translateDna(cdsNoStop);
  const restriction = findRestrictionSites(cdsNoStop);
  const construct = buildConstruct(cdsNoStop, config);
  const parts = annotate(cdsNoStop, config);
  const primers = makePrimers(cdsNoStop, config);
  const mismatch = translated !== proteinState.targetProtein;

  return { ...proteinState, config, cdsNoStop, translated, restriction, construct, parts, primers, mismatch };
}

function setStatus(node, text, type = "neutral") {
  node.textContent = text;
  node.className = `status-pill ${type}`;
}

function showMessage(text) {
  const box = el("messageBox");
  box.hidden = !text;
  box.textContent = text;
}

function renderProteinState() {
  setStatus(el("inputStatus"), `${proteinState.designedProtein.length} aa`, "ok");
  el("proteinDesignSummary").textContent = proteinState.notes.join(" / ") || "No signal peptide or TM edit applied";
  el("reverseTranslateBtn").disabled = false;
  el("copyProteinBtn").disabled = false;
  el("copyDesignedProteinBtn").disabled = false;
  el("finalAnalyzeBtn").disabled = false;
  renderSequenceTab();
}

function renderRestriction(rows) {
  const table = el("restrictionTable");
  table.innerHTML = "";
  const problemRows = rows.filter((row) => row.positions.length);
  setStatus(el("restrictionSummary"), problemRows.length ? `${problemRows.length}개 검출` : "검출 없음", problemRows.length ? "warn" : "ok");
  for (const row of rows) {
    const div = document.createElement("div");
    div.className = "table-row";
    div.innerHTML = `
      <strong>${row.enzyme}</strong>
      <span>${row.site}</span>
      <span class="positions">${row.positions.length ? `위치 ${row.positions.join(", ")}` : "없음"}</span>
    `;
    table.appendChild(div);
  }
}

function renderPrimers(primers) {
  el("primerOutput").className = "primer-output";
  el("primerOutput").innerHTML = `
    <div class="primer-card">
      <h3>${primers.forward.name}</h3>
      <div class="primer-seq">${primers.forward.sequence}</div>
      <div class="detail-text">mode: ${primers.mode} / 5' tail: <code>${primers.forward.tail || "-"}</code> / anneal: <code>${primers.forward.anneal}</code> / Tm ${primers.forward.tm.toFixed(1)} C / GC ${primers.forward.gc.toFixed(1)}%</div>
    </div>
    <div class="primer-card">
      <h3>${primers.reverse.name}</h3>
      <div class="primer-seq">${primers.reverse.sequence}</div>
      <div class="detail-text">mode: ${primers.mode} / 5' tail: <code>${primers.reverse.tail || "-"}</code> / anneal: <code>${primers.reverse.anneal}</code> / Tm ${primers.reverse.tm.toFixed(1)} C / GC ${primers.reverse.gc.toFixed(1)}%</div>
    </div>
  `;
  setStatus(el("primerSummary"), "생성 완료", "ok");
}

function reportText() {
  if (!finalState) return "최종 CDS 분석 후 report가 표시됩니다.";
  const hostLabel = HOST_CONFIG[finalState.config.host].label;
  const problems = finalState.restriction.filter((row) => row.positions.length);
  const vectorSp = vectorEncodedSignalPeptide(finalState.config);
  const strategyLabel = {
    none: "Cytosol",
    native: "Native SP",
    heterologous: "Heterologous SP",
  }[finalState.config.signalStrategy];
  return [
    "Cloning report",
    "==============",
    `Accession no.: ${finalState.config.accession || "-"}`,
    `Purpose: ${finalState.config.purpose || "-"}`,
    `Special notes: ${finalState.config.notes || "-"}`,
    `Host: ${hostLabel}`,
    `Host organism for codon optimization: ${finalState.config.hostOrganism}`,
    `Codon optimization tool: ${finalState.config.codonTool === "thermofisher" ? "Thermo Fisher" : "VectorBuilder"}`,
    `Transfer vector: ${finalState.config.vector}`,
    `Signal peptide strategy: ${strategyLabel}`,
    `Native SP label: ${finalState.config.nativeSignalLabel || "-"}`,
    `SignalP cleavage site probability: ${finalState.config.cleavageProbability || "-"}`,
    `Selected/vector SP: ${vectorSp ? `${vectorSp.label} (${vectorSp.aa}) - vector encoded; GS linker after SP` : finalState.config.signalStrategy === "heterologous" ? `${finalState.config.signalPeptide.label} (${finalState.config.signalPeptide.aa})` : "-"}`,
    `C-terminal tag: ${finalState.config.host === "Baculovirus" ? "vector-encoded C-term His tag; excluded from insert" : `${HIS6_AA} (${HIS6_DNA})`}`,
    `Design notes: ${finalState.notes.join(" / ") || "none"}`,
    `Designed protein length: ${finalState.designedProtein.length} aa`,
    `Final CDS/target DNA length: ${finalState.cdsNoStop.length} bp`,
    `Cloning construct length: ${finalState.construct.length} bp`,
    `Restriction sites: ${problems.length ? problems.map((row) => `${row.enzyme} at ${row.positions.join(",")}`).join("; ") : "none detected in final CDS"}`,
    `Translation check: ${finalState.mismatch ? "warning - optimized DNA translation differs from designed protein" : "OK"}`,
    `Primer mode: ${finalState.primers.mode}`,
  ].join("\n");
}

function renderSequenceTab() {
  const out = el("sequenceOutput");
  if (activeTab === "report") out.textContent = reportText();
  if (activeTab === "protein") out.textContent = proteinState ? fasta("designed_protein", proteinState.designedProtein) : "아미노산 서열 적용 후 표시됩니다.";
  if (activeTab === "cds") out.textContent = finalState ? fasta("final_cds", finalState.cdsNoStop) : "최종 CDS 분석 후 표시됩니다.";
  if (activeTab === "construct") out.textContent = finalState ? fasta("cloning_construct", finalState.construct) : "최종 CDS 분석 후 표시됩니다.";
  if (activeTab === "parts") {
    out.textContent = finalState
      ? finalState.parts.map((part) => `${String(part.start).padStart(5, " ")}-${String(part.end).padStart(5, " ")}  ${part.label.padEnd(16, " ")}  ${part.seq}`).join("\n")
      : "최종 CDS 분석 후 표시됩니다.";
  }
}

function populateVectorOptions() {
  const host = el("hostSelect").value;
  el("vectorSelect").innerHTML = "";
  for (const vector of HOST_CONFIG[host].vectors) {
    const option = document.createElement("option");
    option.value = vector;
    option.textContent = vector;
    el("vectorSelect").appendChild(option);
  }
}

function populateSignalPeptideOptions() {
  const host = el("hostSelect").value;
  el("signalPeptideSelect").innerHTML = "";
  for (const signal of HOST_CONFIG[host].signalPeptides) {
    const option = document.createElement("option");
    option.value = signal.id;
    option.textContent = `${signal.label} (${signal.aa})`;
    el("signalPeptideSelect").appendChild(option);
  }
}

function syncControls() {
  const config = currentConfig();
  const hostLabel = HOST_CONFIG[config.host].label;
  el("signalPeptideSelect").disabled = config.signalStrategy !== "heterologous";
  el("infusionOptions").hidden = !config.infusionEnabled;
  el("forwardOverlap").placeholder = `${config.vector} forward overlap sequence`;
  el("reverseOverlap").placeholder = `${config.vector} reverse overlap sequence`;
  setStatus(el("systemStatus"), `${hostLabel} / ${config.vector}`, "neutral");
}

function syncHostOrganismToHost() {
  const host = el("hostSelect").value;
  if (host === "HEK") el("hostOrganismSelect").value = "Homo sapiens";
  if (host === "CHO") el("hostOrganismSelect").value = "Cricetulus griseus";
  if (host === "Baculovirus") el("hostOrganismSelect").value = "Spodoptera frugiperda";
}

function syncSignalSelectionToVector() {
  const vector = el("vectorSelect").value;
  if (vector.includes("(HM)")) el("signalPeptideSelect").value = "hm";
  if (vector.includes("(SP12)")) el("signalPeptideSelect").value = "sp12";
}

function prepareProtein() {
  const inputProtein = cleanProtein(el("proteinInput").value);
  if (!inputProtein) throw new Error("유효한 아미노산 서열이 없습니다.");
  if (!inputProtein.startsWith("M")) throw new Error("아미노산 서열은 시작 Met(M)을 포함해야 합니다.");
  proteinState = designProtein(inputProtein, currentConfig());
  finalState = null;
  renderProteinState();
  showMessage("");
}

function reverseTranslateCurrent() {
  if (!proteinState) throw new Error("먼저 아미노산 서열을 적용하세요.");
  proteinState.reverseTranslatedDna = reverseTranslateDesignedProtein(proteinState);
  el("reverseTranslationOutput").textContent = fasta("reverse_translated_cds", proteinState.reverseTranslatedDna);
  setStatus(el("reverseStatus"), `${proteinState.reverseTranslatedDna.length} bp`, "ok");
  el("copyReverseTranslatedBtn").disabled = false;
  el("finalAnalyzeBtn").disabled = false;
}

function renderFinal() {
  finalState = finalAnalyze();
  renderRestriction(finalState.restriction);
  renderPrimers(finalState.primers);
  el("copyForwardPrimerBtn").disabled = false;
  el("copyReversePrimerBtn").disabled = false;
  el("copyPrimerPairBtn").disabled = false;
  if (finalState.mismatch) {
    showMessage("경고: 입력한 codon optimization DNA의 번역 결과가 설계 단백질과 다릅니다.");
  } else {
    showMessage("");
  }
  renderSequenceTab();
}

async function copyText(text) {
  await navigator.clipboard.writeText(text);
}

function resetApp() {
  proteinState = null;
  finalState = null;
  el("proteinInput").value = "";
  el("optimizedDnaInput").value = "";
  el("reverseTranslationOutput").textContent = "DNA 서열은 한 줄에 105 bp씩 표시됩니다.";
  el("restrictionTable").textContent = "최종 DNA 분석 후 표시됩니다.";
  el("restrictionTable").className = "table-like empty-state";
  el("primerOutput").textContent = "최종 CDS 분석 후 권장 primer가 표시됩니다.";
  el("primerOutput").className = "primer-output empty-state";
  setStatus(el("inputStatus"), "대기", "neutral");
  setStatus(el("reverseStatus"), "대기", "neutral");
  setStatus(el("restrictionSummary"), "대기", "neutral");
  setStatus(el("primerSummary"), "대기", "neutral");
  ["copyProteinBtn", "copyDesignedProteinBtn", "reverseTranslateBtn", "copyReverseTranslatedBtn", "finalAnalyzeBtn", "copyForwardPrimerBtn", "copyReversePrimerBtn", "copyPrimerPairBtn"].forEach((id) => {
    el(id).disabled = true;
  });
  el("proteinDesignSummary").textContent = "아미노산 서열 적용 후 signal peptide/TM 처리 결과가 표시됩니다.";
  showMessage("");
  renderSequenceTab();
}

document.addEventListener("DOMContentLoaded", () => {
  populateVectorOptions();
  populateSignalPeptideOptions();
  syncHostOrganismToHost();
  syncControls();

  el("hostSelect").addEventListener("change", () => {
    populateVectorOptions();
    populateSignalPeptideOptions();
    syncHostOrganismToHost();
    syncSignalSelectionToVector();
    syncControls();
  });
  el("vectorSelect").addEventListener("change", () => {
    syncSignalSelectionToVector();
    syncControls();
  });
  ["vectorSelect", "signalStrategySelect", "signalPeptideSelect", "nativeSignalAa", "cleavageProbability", "nativeSignalLabel", "tmStart", "tmEnd", "infusionToggle", "forwardOverlap", "reverseOverlap", "codonToolSelect", "hostOrganismSelect", "accessionInput", "purposeInput", "notesInput"].forEach((id) => {
    el(id).addEventListener("input", syncControls);
    el(id).addEventListener("change", syncControls);
  });

  el("prepareProteinBtn").addEventListener("click", () => {
    try { prepareProtein(); } catch (error) { showMessage(error.message); setStatus(el("inputStatus"), "오류", "danger"); }
  });
  el("previewProteinDesignBtn").addEventListener("click", () => {
    try { prepareProtein(); } catch (error) { showMessage(error.message); }
  });
  el("applyProteinDesignBtn").addEventListener("click", () => {
    try { prepareProtein(); } catch (error) { showMessage(error.message); }
  });
  el("reverseTranslateBtn").addEventListener("click", () => {
    try { reverseTranslateCurrent(); } catch (error) { showMessage(error.message); }
  });
  el("finalAnalyzeBtn").addEventListener("click", () => {
    try { renderFinal(); } catch (error) { showMessage(error.message); }
  });
  el("optimizedDnaInput").addEventListener("input", () => {
    if (proteinState) el("finalAnalyzeBtn").disabled = false;
  });

  el("copyProteinBtn").addEventListener("click", () => proteinState && copyText(fasta("designed_protein", proteinState.designedProtein)));
  el("copyDesignedProteinBtn").addEventListener("click", () => proteinState && copyText(fasta("designed_protein", proteinState.designedProtein)));
  el("copyReverseTranslatedBtn").addEventListener("click", () => proteinState?.reverseTranslatedDna && copyText(formatSeq(proteinState.reverseTranslatedDna)));
  el("copyForwardPrimerBtn").addEventListener("click", () => finalState && copyText(`${finalState.primers.forward.name}\n${finalState.primers.forward.sequence}`));
  el("copyReversePrimerBtn").addEventListener("click", () => finalState && copyText(`${finalState.primers.reverse.name}\n${finalState.primers.reverse.sequence}`));
  el("copyPrimerPairBtn").addEventListener("click", () => finalState && copyText([
    `${finalState.primers.forward.name}: ${finalState.primers.forward.sequence}`,
    `${finalState.primers.reverse.name}: ${finalState.primers.reverse.sequence}`,
  ].join("\n")));

  el("loadExampleBtn").addEventListener("click", () => {
    el("proteinInput").value = ">example_protein\nMKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAE";
    el("nativeSignalAa").value = "MKWVTFISLLFLFSSAYS";
    el("nativeSignalLabel").value = "native HSA-like SP 1-18 aa";
    prepareProtein();
  });
  el("clearBtn").addEventListener("click", resetApp);

  document.querySelectorAll(".tab").forEach((button) => {
    button.addEventListener("click", () => {
      document.querySelectorAll(".tab").forEach((tab) => tab.classList.remove("active"));
      button.classList.add("active");
      activeTab = button.dataset.tab;
      renderSequenceTab();
    });
  });
});
