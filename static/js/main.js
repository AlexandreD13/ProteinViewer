const canvas = document.getElementById("glCanvas");
const gl = canvas.getContext("webgl");

if (!gl) {
    alert("WebGL not supported!");
}

let showAtoms = true;
let showBonds = true;
let atomSize = 25.0;

function resizeCanvas() {
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
    gl.viewport(0, 0, canvas.width, canvas.height);
}

resizeCanvas();
window.addEventListener("resize", resizeCanvas);

fetch("/protein/6lcw.pdb")
    .then(response => response.text())
    .then(parsePDB)
    .then(initWebGL);

function parsePDB(pdbData) {
    let atomPositions = [];
    let atomicNumbers = [];
    let bonds = [];

    const atomMap = {"H": 1.0, "C": 6.0, "N": 7.0, "O": 8.0, "S": 16.0};

    let lines = pdbData.split("\n");
    let atoms = [];

    for (let line of lines) {
        if (line.startsWith("ATOM")) {
            let x = parseFloat(line.substring(30, 38));
            let y = parseFloat(line.substring(38, 46));
            let z = parseFloat(line.substring(46, 54));

            let atomSymbol = line.substring(76, 78).trim();
            let atomicNumber = atomMap[atomSymbol] || 0.0;

            atoms.push({x, y, z, index: atoms.length});
            atomPositions.push(x, y, z);
            atomicNumbers.push(atomicNumber);
        }
    }

    // Bond detection based on distance threshold
    const bondThreshold = 2.0; // Approximate covalent bond length
    for (let i = 0; i < atoms.length; i++) {
        for (let j = i + 1; j < atoms.length; j++) {
            let dx = atoms[i].x - atoms[j].x;
            let dy = atoms[i].y - atoms[j].y;
            let dz = atoms[i].z - atoms[j].z;
            let distance = Math.sqrt(dx * dx + dy * dy + dz * dz);

            if (distance < bondThreshold) {
                bonds.push(i, j);
            }
        }
    }

    // Update the protein info on the page
    document.getElementById("proteinName").textContent = "Example Protein Name"; // You can fetch this from PDB data if needed
    document.getElementById("atomCount").textContent = atoms.length;
    document.getElementById("bondCount").textContent = bonds.length;

    return {
        positions: new Float32Array(atomPositions),
        atomicNumbers: new Float32Array(atomicNumbers),
        bonds: new Uint16Array(bonds),
    };
}

function initWebGL(atomData) {
    let rotationX = 0, rotationY = 0, zoom = 2.0;
    let translationX = 0, translationY = 0;

    const projectionMatrix = mat4.create();
    const modelViewMatrix = mat4.create();

    const vsSource = `
                attribute vec3 atomPosition;
                attribute float atomicNumber;
                uniform mat4 uModelViewMatrix;
                uniform mat4 uProjectionMatrix;
                uniform float uZoom;
                uniform float uAtomSize;
                varying float vAtomicNumber;

                void main() {
                    vAtomicNumber = atomicNumber;
                    gl_PointSize = uAtomSize / uZoom;
                    gl_Position = uProjectionMatrix * uModelViewMatrix * vec4(atomPosition.xyz / 50.0, 1.0);
                }
            `;

    const fsSource = `
                precision mediump float;
                varying float vAtomicNumber;

                vec4 getColor(float atomicNumber) {
                    if (abs(atomicNumber - 1.0) < 0.01) return vec4(1.0, 1.0, 1.0, 1.0);    // Hydrogen (White)
                    if (abs(atomicNumber - 6.0) < 0.01) return vec4(1.0, 0.0, 0.0, 1.0);    // Carbon (Red)
                    if (abs(atomicNumber - 7.0) < 0.01) return vec4(0.0, 0.0, 1.0, 1.0);    // Nitrogen (Blue)
                    if (abs(atomicNumber - 8.0) < 0.01) return vec4(0.0, 1.0, 0.0, 1.0);    // Oxygen (Green)
                    if (abs(atomicNumber - 16.0) < 0.01) return vec4(1.0, 1.0, 0.0, 1.0);   // Sulfur (Yellow)
                    return vec4(0.5, 0.5, 0.5, 1.0);                                        // Default gray color
                }

                void main() {
                    vec2 coord = gl_PointCoord - vec2(0.5);
                    float r2 = dot(coord, coord);
                    if (r2 > 0.25) discard;

                    gl_FragColor = getColor(vAtomicNumber);
                }
            `;

    const bondVsSource = `
                attribute vec3 bondPosition;
                uniform mat4 uModelViewMatrix;
                uniform mat4 uProjectionMatrix;

                void main() {
                    gl_Position = uProjectionMatrix * uModelViewMatrix * vec4(bondPosition.xyz / 50.0, 1.0);
                }
            `;

    const bondFsSource = `
                precision mediump float;

                void main() {
                    gl_FragColor = vec4(0.8, 0.8, 0.8, 1.0);
                }
            `;

    const vertexShader = gl.createShader(gl.VERTEX_SHADER);
    gl.shaderSource(vertexShader, vsSource);
    gl.compileShader(vertexShader);

    const fragmentShader = gl.createShader(gl.FRAGMENT_SHADER);
    gl.shaderSource(fragmentShader, fsSource);
    gl.compileShader(fragmentShader);

    const shaderProgram = gl.createProgram();
    gl.attachShader(shaderProgram, vertexShader);
    gl.attachShader(shaderProgram, fragmentShader);
    gl.linkProgram(shaderProgram);
    gl.useProgram(shaderProgram);

    const bondVertexShader = gl.createShader(gl.VERTEX_SHADER);
    gl.shaderSource(bondVertexShader, bondVsSource);
    gl.compileShader(bondVertexShader);

    const bondFragmentShader = gl.createShader(gl.FRAGMENT_SHADER);
    gl.shaderSource(bondFragmentShader, bondFsSource);
    gl.compileShader(bondFragmentShader);

    const bondShaderProgram = gl.createProgram();
    gl.attachShader(bondShaderProgram, bondVertexShader);
    gl.attachShader(bondShaderProgram, bondFragmentShader);
    gl.linkProgram(bondShaderProgram);

    const positionAttrib = gl.getAttribLocation(shaderProgram, "atomPosition");
    const atomicNumberAttrib = gl.getAttribLocation(shaderProgram, "atomicNumber");
    gl.enableVertexAttribArray(positionAttrib);
    gl.enableVertexAttribArray(atomicNumberAttrib);

    const uModelViewMatrix = gl.getUniformLocation(shaderProgram, "uModelViewMatrix");
    const uProjectionMatrix = gl.getUniformLocation(shaderProgram, "uProjectionMatrix");
    const zoomUniform = gl.getUniformLocation(shaderProgram, "uZoom");
    const atomSizeUniform = gl.getUniformLocation(shaderProgram, "uAtomSize");

    mat4.perspective(projectionMatrix, Math.PI / 4, canvas.width / canvas.height, 0.1, 100);
    gl.uniformMatrix4fv(uProjectionMatrix, false, projectionMatrix);

    const positionBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, atomData.positions, gl.STATIC_DRAW);
    gl.vertexAttribPointer(positionAttrib, 3, gl.FLOAT, false, 0, 0);

    const bondBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, bondBuffer);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, atomData.bonds, gl.STATIC_DRAW);

    const atomicNumberBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, atomicNumberBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, atomData.atomicNumbers, gl.STATIC_DRAW);
    gl.vertexAttribPointer(atomicNumberAttrib, 1, gl.FLOAT, false, 0, 0);

    // Mouse and Zoom Controls
    let isDragging = false;
    let lastX, lastY;
    let isCtrlPressed = false;

    canvas.addEventListener("mousedown", (event) => {
        isDragging = true;
        lastX = event.clientX;
        lastY = event.clientY;
        isCtrlPressed = event.ctrlKey;
    });

    canvas.addEventListener("mousemove", (event) => {
        if (isDragging) {
            let deltaX = event.clientX - lastX;
            let deltaY = event.clientY - lastY;

            if (isCtrlPressed) {
                translationX += deltaX * 0.01;
                translationY -= deltaY * 0.01;
            } else {
                rotationY += deltaX * 0.01;
                rotationX += deltaY * 0.01;
            }

            lastX = event.clientX;
            lastY = event.clientY;

            updateMatrix();
        }
    });

    canvas.addEventListener("mouseup", () => {
        isDragging = false;
        isCtrlPressed = false;
    });

    canvas.addEventListener("mouseleave", () => {
        isDragging = false;
        isCtrlPressed = false;
    });

    canvas.addEventListener("wheel", (event) => {
        zoom += event.deltaY * 0.01;
        zoom = Math.max(0.1, Math.min(zoom, 20.0));
        updateMatrix();
    });

    function updateMatrix() {
        mat4.identity(modelViewMatrix);
        mat4.translate(modelViewMatrix, modelViewMatrix, [translationX, translationY, -zoom]);
        mat4.rotateX(modelViewMatrix, modelViewMatrix, rotationX);
        mat4.rotateY(modelViewMatrix, modelViewMatrix, rotationY);

        gl.useProgram(shaderProgram);
        gl.uniformMatrix4fv(uModelViewMatrix, false, modelViewMatrix);
        gl.uniformMatrix4fv(uProjectionMatrix, false, projectionMatrix);
        gl.uniform1f(zoomUniform, zoom);
        gl.uniform1f(atomSizeUniform, atomSize);
    }

    function render() {
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
        updateMatrix();

        // Draw atoms only if showAtoms is true
        if (showAtoms) {
            gl.useProgram(shaderProgram);
            gl.uniformMatrix4fv(uModelViewMatrix, false, modelViewMatrix);
            gl.uniformMatrix4fv(uProjectionMatrix, false, projectionMatrix);
            gl.drawArrays(gl.POINTS, 0, atomData.positions.length / 3);
        }

        // Draw bonds only if showBonds is true
        if (showBonds) {
            gl.useProgram(bondShaderProgram);
            gl.uniformMatrix4fv(gl.getUniformLocation(bondShaderProgram, "uModelViewMatrix"), false, modelViewMatrix);
            gl.uniformMatrix4fv(gl.getUniformLocation(bondShaderProgram, "uProjectionMatrix"), false, projectionMatrix);
            gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, bondBuffer);
            gl.drawElements(gl.LINES, atomData.bonds.length, gl.UNSIGNED_SHORT, 0);
        }

        requestAnimationFrame(render);
    }

    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.enable(gl.DEPTH_TEST);
    render();

    document.getElementById("toggleAtoms").addEventListener("click", () => {
        showAtoms = !showAtoms;
        document.getElementById("toggleAtoms").textContent = showAtoms ? "Hide Atoms" : "Show Atoms";
    });

    document.getElementById("toggleBonds").addEventListener("click", () => {
        showBonds = !showBonds;
        document.getElementById("toggleBonds").textContent = showBonds ? "Hide Bonds" : "Show Bonds";
    });

    document.getElementById("atomSizeSlider").addEventListener("input", (event) => {
        atomSize = parseFloat(event.target.value);
        document.getElementById("atomSizeValue").textContent = atomSize;
        updateMatrix();
    });
}