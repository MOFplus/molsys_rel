// parse input
var edges = process.argv[2].split("|");
var path = process.argv[3];

// construct lqg in systreKey format (array of edges [[v1 v2 [1 0 0]] ...
var lqg = [];
for (let i = 0; i < edges.length; i++) {
    let numbers = edges[i].split(" ");
    numbers = numbers.map(Number);
    let rel = numbers.slice(2,5);
    lqg.push([numbers[0], numbers[1], rel]);
}

// load package, assuming that systreKey lives in the same dir as molsys (usually the git repository dir)
const sk = require(path + '/dist/systreKey').systreKey;

// get key and write it to stdout
let out = sk(lqg);
console.log(out);



