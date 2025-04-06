const express = require('express');
const cors = require('cors');
const path = require('path');
const { parsePDBFile } = require('./parser');
require('dotenv').config();

const app = express();
const PORT = process.env.PORT || 5000;

app.use(cors());
app.use(express.json());

app.get('/api/protein/:filename', (req, res) => {
	const filename = req.params.filename;

	try {
		const parsedData = parsePDBFile(filename);
		res.json(parsedData);
	} catch (err) {
		console.error(err);
		res.status(404).json({ error: 'File not found or could not be parsed.' });
	}
});

app.listen(PORT, () => {
	console.log(`Server is running on http://localhost:${PORT}`);
});
