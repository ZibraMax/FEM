function readSingleFile(e) {
    var file = e.target.files[0];
    console.log(file, e);
    if (!file) {
        return;
    }
    var reader = new FileReader();
    reader.onload = function (e) {
        var contents = e.target.result;
        loadMSH(contents);
    };
    reader.readAsText(file);
}

async function loadMSH(str) {
    let lineas = str.split("\n");
    let parametros = lineas[0].split("\t").map((x) => parseFloat(x));
    for (let i = parametros.length - 1; i < 9; i++) {
        parametros.push(0);
    }
    let param = {
        ngdl: parametros[0],
        nel: parametros[1],
        nseg: parametros[2],
        ncbe: parametros[3],
        ncbn: parametros[4],
        nvn: parametros[5],
        nholes: parametros[6],
        nfillets: parametros[7],
        nmask: parametros[8],
    };
    let gdls = [];
    let elements = [];
    let types = [];
    let segments = [];
    let cbe = [];
    let cbn = [];
    let holes = [];
    let fillets = [];
    let mask = [];
    let maxx = 0;
    let maxy = 0;
    let j = 1;
    let cx = 0;
    let cy = 0;
    for (let i = 0; i < param.ngdl; i++) {
        let coord = lineas[i + j].split("\t").map((x) => parseFloat(x));
        gdls.push(coord);
        let x = coord[0];
        cx += x;
        let y = coord[1];
        cy += y;
        if (x > maxx) {
            maxx = x;
        }
        if (y > maxy) {
            maxy = y;
        }
    }
    cx = cx / param.ngdl;
    cy = cy / param.ngdl;
    j += param.ngdl;
    for (let i = 0; i < param.nel; i++) {
        types.push(lineas[i + j]);
    }
    j += param.nel;
    for (let i = 0; i < param.nel; i++) {
        elements.push(lineas[i + j].split("\t").map((x) => parseInt(x)));
    }
    j += param.nel;
    for (let i = 0; i < param.nseg; i++) {
        segments.push(lineas[i + j].split("\t").map((x) => parseInt(x)));
    }
    j += param.nseg;
    for (let i = 0; i < param.ncbe; i++) {
        cbe.push(lineas[i + j].split("\t").map((x) => parseFloat(x)));
    }
    j += param.ncbe;
    for (let i = 0; i < param.ncbn; i++) {
        cbn.push(lineas[i + j].split("\t").map((x) => parseFloat(x)));
    }
    j += param.ncbn;

    for (let i = 0; i < param.nholes; i++) {
        holes.push(eval("(" + lineas[i + j] + ")"));
    }
    j += param.nholes;
    for (let i = 0; i < param.nfillets; i++) {
        fillets.push(eval("(" + lineas[i + j] + ")"));
    }
    j += param.nfillets;
    for (let i = 0; i < param.nmask; i++) {
        mask.push(lineas[i + j].split("\t").map((x) => parseFloat(x)));
    }
    j += param.nmask;

    let geometry = {
        gdls: gdls,
        elements: elements,
        types: types,
        segments: segments,
        cbe: cbe,
        cbn: cbn,
        holes: holes,
        fillets: fillets,
        mask: mask,
    };
    GEOMETRY = geometry;
    let multx = two.width / maxx;
    let multy = two.height / maxy;
    let mult = Math.min(multx, multy);
    console.log(mult, cx, cy);

    draw(GEOMETRY, mult, cx, cy);
    two.unbind("update");
}

function drawCircle(coord, mult) {
    return new Promise((resolve) => {
        const circle = two.makeCircle(
            coord[0] * mult,
            two.height - coord[1] * mult,
            3
        );
        circle.fill = "orangered";
        circle.noStroke();
        resolve(circle);
    });
}
function drawLine(coord, coordi, mult) {
    return new Promise((resolve) => {
        const line = two.makeLine(
            coord[0] * mult,
            two.height - coord[1] * mult,
            coordi[0] * mult,
            two.height - coordi[1] * mult
        );
        line.stroke = "#f7924a";
        line.linewidth = 2;
        resolve(line);
    });
}

let nn = [];
let nnl = [];
async function draw(geometry, mult, cx, cy) {
    var lineas = two.makeGroup(nnl);
    var nodos = two.makeGroup(nn);
    console.time("Contando");
    for (let i = 0; i < geometry.elements.length; i++) {
        for (let j = 0; j < 3; j++) {
            let coord = geometry.gdls[geometry.elements[i][j]];
            let coordi = geometry.gdls[geometry.elements[i][j + 1]];
            if (j == 2) {
                coordi = geometry.gdls[geometry.elements[i][0]];
            }
            let circle = await drawCircle(coord, mult);
            let line = await drawLine(coord, coordi, mult);
            nodos.add(circle);
            lineas.add(line);
        }
    }
    let scale = 0.75;
    lineas.scale = scale;
    lineas.translation.set(
        two.width / 2 - cx * mult * scale,
        two.height / 2 - cy * mult * scale
    );
    nodos.scale = scale;
    nodos.translation.set(
        two.width / 2 - cx * mult * scale,
        two.height / 2 - cy * mult * scale
    );
    two.update();
    console.timeEnd("Contando");
}

var GEOMETRY = undefined;
var elem = document.getElementById("main");
let pad = 2 * parseInt(getComputedStyle(elem).padding, 10);
let width = elem.offsetWidth - pad;
let height = elem.offsetHeight - pad;
var two = new Two({
    autostart: true,
    width: width,
    height: height,
    type: Two.Types.canvas,
}).appendTo(elem);
var mouse = new Two.Vector();
two.renderer.domElement.addEventListener("mousedown", mousedown, false);
two.renderer.domElement.addEventListener("mousewheel", mousewheel, false);
var zui = new Two.ZUI(two.scene);
zui.addLimits(0.06, 8);
var $stage = $(two.renderer.domElement);
var SELECCIONANDO = false;

function toogleSelect() {
    SELECCIONANDO = !SELECCIONANDO;
}

function mousedown(e) {
    mouse.x = e.clientX;
    mouse.y = e.clientY;
    window.addEventListener("mousemove", mousemove, false);
    window.addEventListener("mouseup", mouseup, false);
}

function mousemove(e) {
    var dx = e.clientX - mouse.x;
    var dy = e.clientY - mouse.y;
    zui.translateSurface(dx, dy);
    mouse.set(e.clientX, e.clientY);
}

function mouseup(e) {
    window.removeEventListener("mousemove", mousemove, false);
    window.removeEventListener("mouseup", mouseup, false);
}

function mousewheel(e) {
    var dy = (e.wheelDeltaY || -e.deltaY) / 1000;
    zui.zoomBy(dy, e.clientX, e.clientY);
}

function panstart(e) {
    var touch = e.touches[0];
    mouse.x = touch.clientX;
    mouse.y = touch.clientY;
}

function panmove(e) {
    var touch = e.touches[0];
    var dx = touch.clientX - mouse.x;
    var dy = touch.clientY - mouse.y;
    zui.translateSurface(dx, dy);
    mouse.set(touch.clientX, touch.clientY);
}

$stage.bind("mousewheel wheel", function (event) {
    var e = event.originalEvent;
    e.stopPropagation();
    e.preventDefault();

    var dy = (e.wheelDeltaY || -e.deltaY) / 1000;

    zui.zoomBy(dy, e.clientX, e.clientY);
});
document
    .getElementById("file-input")
    .addEventListener("change", readSingleFile, false);
