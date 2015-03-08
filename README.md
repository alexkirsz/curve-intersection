# curve-intersection

Bezier curve intersection algorithm and utilities.

Extracted from the [Paper.js](https://github.com/paperjs/paper.js) implementation of bezier clipping.

## Installation

```
npm install --save curve-intersection
```

## Usage

If your platform doesn't support es6 yet, you can `require('curve-intersection/es3')`.

```javascript
import {curveIntersections} from 'curve-intersection';

// coordinates of the control points
let curves = [
  [ 25.3, 21.4, -93.4, -180.5, 90.9, 177.2, -31, -15.8 ],
  [ 26.9, -22.6, -196.3, 48.300000000000004, 193.4, -52, -21.8, 24 ]
];

let intersections = curveIntersections(curves[0], curves[1]);

let canvas = document.createElement('canvas');
canvas.width = 200;
canvas.height = 200;
document.body.appendChild(canvas);
let ctx = canvas.getContext('2d');
ctx.translate(100, 100);

for (let curve of curves) {
  ctx.beginPath();
  ctx.moveTo(curve[0], curve[1]);
  ctx.bezierCurveTo(curve[2], curve[3], curve[4], curve[5], curve[6], curve[7]);
  ctx.closePath();
  ctx.stroke();
}

for (let intersection of intersections) {
  let [t1, p1, t2, p2] = intersection;
  ctx.beginPath();
  ctx.arc(p1[0], p1[1], 5, 0, 2 * Math.PI, false);
  ctx.closePath();
  ctx.fill();
}
```
