const express = require('express');
const cors = require('cors');
const bodyParser = require('body-parser');

const { parsePDBFile } = require('./parser');

require('dotenv').config();

const app = express();
const PORT = process.env.PORT || 5000;

app.use(bodyParser.json({ limit: '10mb' }));
app.use(bodyParser.urlencoded({ limit: '10mb', extended: true }));
app.use(cors());
app.use(express.json());

app.post('/api/protein/upload', (req, res) => {
	const { name, data } = req.body;

	if (!data) {
		return res.status(400).json({ error: 'No file content provided.' });
	}

	try {
		const parsed = parsePDBFile(name, data);
		res.json(parsed);
	} catch (err) {
		console.error(err);
		res.status(500).json({ error: 'Failed to parse PDB content.' });
	}
});

app.listen(PORT, () => {
	console.log(`Server is running on http://localhost:${PORT}`);
});
