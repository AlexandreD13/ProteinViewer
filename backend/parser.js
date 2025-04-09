const fs = require("fs");
const path = require("path");

const atomMap = { H: 1.0, C: 6.0, N: 7.0, O: 8.0, S: 16.0 };

const residueMap = {
	GLY: "GLYCINE", ALA: "ALANINE", VAL: "VALINE", LEU: "LEUCINE",
	ILE: "ISOLEUCINE", SER: "SERINE", THR: "THREONINE", TYR: "TYROSINE",
	CYS: "CYSTEINE", MET: "METHIONINE", ASP: "ASPARTIC ACID", ASN: "ASPARAGINE",
	GLU: "GLUTAMIC ACID", GLN: "GLUTAMINE", ARG: "ARGININE", LYS: "LYSINE",
	HIS: "HISTIDINE", PHE: "PHENYLALANINE", TRP: "TRYPTOPHAN", PRO: "PROLINE"
};

function parseRecord(recordsString) {
	return recordsString
		.split("MOL_ID:")
		.slice(1)
		.reduce((acc, entry) => {
			const parts = entry.trim().split(";");
			const compoundId = parts[0].trim();
			const record = {};

			parts.forEach(part => {
				const [key, value] = part.split(":").map(s => s.trim());
				if (key && value) record[key] = value;
			});

			if (compoundId) acc[compoundId] = record;
			return acc;
		}, {});
}

function parseBONDS(atoms) {
	const bonds = [];
	const thresholdSquared = 2.0 * 2.0;

	for (let i = 0; i < atoms.length - 1; i++) {
		const a = atoms[i];
		for (let j = i + 1; j < atoms.length; j++) {
			const b = atoms[j];
			const dx = a.x - b.x;
			const dy = a.y - b.y;
			const dz = a.z - b.z;
			const distSquared = dx * dx + dy * dy + dz * dz;

			if (distSquared < thresholdSquared) bonds.push(i, j);
		}
	}

	return bonds;
}

function parseREVDAT(line, revData) {
	const modNumber = line.substring(7, 10).trim();
	const modDate = new Date(line.substring(13, 22).trim());
	const modID = line.substring(23, 27).trim();
	const modType = line.at(31);
	const modRecord = line.substring(39).trim().split(" ").filter(s => s !== "");

	revData.push({
		modNumber, modDate, modID, modType, modRecord
	});
}

function parseDBREF(line, dbRef) {
	const id = line.substring(7, 11).trim();
	const chainId = line.at(12);
	const seqBegin = line.substring(14, 18).trim();
	const insertBegin = line.at(18);
	const seqEnd = line.substring(20, 24).trim();
	const insertEnd = line.at(24);
	const database = line.substring(26, 32).trim();
	const dbAccession = line.substring(33, 41).trim();
	const dbIdCode = line.substring(42, 54).trim();
	const dbSeqBegin = line.substring(55, 60).trim();
	const idBNSBeg = line.at(60);
	const dbSeqEnd = line.substring(62, 67).trim();
	const dbINSEnd = line.at(67);

	dbRef.push({
		id, chainId, seqBegin, insertBegin, seqEnd, insertEnd,
		database, dbAccession, dbIdCode, dbSeqBegin, idBNSBeg,
		dbSeqEnd, dbINSEnd
	});
}

function parseSEQADV(line, seqAdv) {
	const idCode = line.substring(7, 11).trim();
	const resName = line.substring(12, 15).trim();
	const chainId = line.at(16);
	const seqNum = line.substring(18, 22).trim();
	const iCode = line.at(22);
	const database = line.substring(24, 28).trim();
	const dbAccession = line.substring(29, 38).trim();
	const dbRes = line.substring(39, 42).trim();
	const dbSeq = line.substring(43, 48).trim();
	const conflict = line.substring(49, 70).trim();

	seqAdv.push({
		idCode, resName, chainId, seqNum, iCode,
		database, dbAccession, dbRes, dbSeq, conflict
	});
}

function parseSEQRES(lines) {
	const chains = {};

	for (const line of lines) {
		if (!line.startsWith('SEQRES')) continue;

		const chainID = line.substring(11, 12);
		const residues = line.substring(19).trim().split(/\s+/);

		if (!chains[chainID]) {
			chains[chainID] = [];
		}

		chains[chainID].push(...residues);
	}

	return chains;
}

function parseHET(line, het) {
	const hetID = line.substring(7, 10).trim();
	const chainID = line.at(12);
	const seqNum = line.substring(13, 17).trim();
	const iCode = line.at(17);
	const numHetAtoms = line.substring(20, 25).trim();
	const text = line.substring(30, 70).trim();

	het.push({
		hetID, chainID, seqNum,
		iCode, numHetAtoms, text
	});
}

function parseATOMS(line, atoms) {
	const name = line.substring(12, 16).trim();
	const residueName = line.substring(17, 20).trim();
	const fullResidueName = residueMap[residueName];
	const chainID = line.substring(21, 22).trim();
	const x = parseFloat(line.substring(30, 38));
	const y = parseFloat(line.substring(38, 46));
	const z = parseFloat(line.substring(46, 54));
	const occupancy = parseFloat(line.substring(54, 60));
	const tempFactor = parseFloat(line.substring(60, 66));

	const element = name[0];
	const atomicNumber = atomMap[element] || 0.0;

	atoms.push({
		name, residueName, fullResidueName, chainID,
		x, y, z, occupancy, tempFactor, atomicNumber
	});
}

function parseJRNL(lines) {
	const journal = {
		authors: [],
		title: "",
		edit: "",
		reference: "",
		publ: "",
		refn: "",
		pmid: "",
		doi: ""
	};

	let currentField = "";

	for (const line of lines) {
		if (!line.startsWith("JRNL")) continue;

		const field = line.substring(12, 16).trim();
		const content = line.substring(19).trim();

		switch (field) {
			case "AUTH":
				journal.authors.push(...content.split(",").map(a => a.trim()));
				currentField = "AUTH";
				break;
			case "TITL":
				journal.title += (journal.title ? " " : "") + content;
				currentField = "TITL";
				break;
			case "EDIT":
				journal.edit = content;
				currentField = "EDIT";
				break;
			case "REF":
				journal.reference = {
					pubName: content.substring(0, 28).trim(),
					volume: content.substring(32, 36).trim(),
					page: content.substring(37, 42).trim(),
					year: content.substring(43, 47).trim()
				};
				currentField = "REF";
				break;
			case "PUBL":
				journal.publ = content;
				currentField = "PUBL";
				break;
			case "REFN":
				journal.refn = content;
				currentField = "REFN";
				break;
			case "PMID":
				journal.pmid = content;
				currentField = "PMID";
				break;
			case "DOI":
				journal.doi = content;
				currentField = "DOI";
				break;
			default:
				if (currentField === "TITL") {
					journal.title += " " + content;
				} else if (currentField === "REF") {
					journal.reference += " " + content;
				}
				break;
		}
	}

	journal.authors = journal.authors.filter(s => s !== "");

	return journal;
}

function parsePDBFile(filename) {
	const filePath = path.join(__dirname, "..", "data", "pdb_files", filename);
	const data = fs.readFileSync(filePath, "utf8");

	let proteinClassification = "";
	let depositedDate = ""
	let proteinName = "";
	let proteinDescription = "";
	let compoundsString = "";
	let sourcesString = "";
	let keywordsString = "";
	let expDataString = "";
	let authorsString = "";

	const atoms = [];
	const revData = [];
	const dbRef = [];
	const seqAdv = [];
	const het = [];

	const lines = data.split("\n");

	for (let line of lines) {
		if (line.startsWith("HEADER")) {
			const filteredLine = line.replace("\r", "").split(" ").filter(Boolean);

			proteinClassification = filteredLine.slice(1, -2).join(" ");
			depositedDate = filteredLine.at(-2);
			proteinName = filteredLine.at(-1);
		}

		if (line.startsWith("TITLE")) {
			proteinDescription += line.substring(10, 80).trimEnd();
		}

		if (line.startsWith("COMPND")) {
			compoundsString += line.substring(10, 80).trimEnd();
		}

		if (line.startsWith("SOURCE")) {
			sourcesString += line.substring(10, 80).trimEnd();
		}

		if (line.startsWith("KEYWDS")) {
			keywordsString += line.substring(10, 80).trimEnd();
		}

		if (line.startsWith("EXPDTA")) {
			expDataString += line.substring(10, 80).trimEnd();
		}

		if (line.startsWith("AUTHOR")) {
			authorsString += line.substring(10, 80).trimEnd();
		}

		if (line.startsWith("REVDAT")) {
			parseREVDAT(line, revData);
		}

		if (line.startsWith("DBREF")) {
			parseDBREF(line, dbRef);
		}

		if (line.startsWith("SEQADV")) {
			parseSEQADV(line, seqAdv);
		}

		if (line.startsWith("HET ")) {
			parseHET(line, het);
		}

		// ^^^^^^^^^^^^^^^^^^

		// parseRemarks();

		// vvvvvvvvvvvvvvvvvv

		if (line.startsWith("ATOM")) {
			parseATOMS(line, atoms);
		}
	}

	const bonds = parseBONDS(atoms);

	return {
		proteinName: proteinName,
		proteinClassification: proteinClassification,
		depositedDate: new Date(depositedDate),
		proteinDescription: proteinDescription,
		compounds: parseRecord(compoundsString),
		sources: parseRecord(sourcesString),
		keywords: keywordsString.split(", "),
		expData: expDataString.split(", "),
		authors: authorsString.split(","),
		revData: revData,
		journal: parseJRNL(lines),
		dbRef: dbRef,
		seqAdv: seqAdv,
		seqres: parseSEQRES(lines),
		het: het,
		atomCount: atoms.length,
		bondCount: bonds.length / 2,
		atoms: atoms,
		bonds: bonds,
	};
}

module.exports = {parsePDBFile};
