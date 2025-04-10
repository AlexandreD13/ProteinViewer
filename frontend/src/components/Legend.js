import React from "react";
import "./Legend.css";

const atomColors = [
	{ element: "H", name: "Hydrogen", color: "#FFFFFF" },
	{ element: "C", name: "Carbon", color: "#FF0000" },
	{ element: "N", name: "Nitrogen", color: "#0000FF" },
	{ element: "O", name: "Oxygen", color: "#00FF00" },
	{ element: "S", name: "Sulfur", color: "#FFFF00" },
];

export default function Legend() {
	return (
		<div id="legend">
			<strong><u>Atom Legend:</u></strong><br/><br/>
			{atomColors.map(({ element, name, color }) => (
				<div key={element} className="legend-item">
					<div
						className="legend-swatch"
						style={{ backgroundColor: color, borderRadius: "50px" }}
					></div>
					<div className="legend-text">
						<strong>{element}</strong> - {name}
					</div>
				</div>
			))}
		</div>
	);
}
