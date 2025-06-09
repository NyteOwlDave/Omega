
// MDN Playground
// https://developer.mozilla.org/en-US/play

// Omega Graphics API
// http://dave-ryzen/omega/omega.html

// BOB Image Viewer
// http://dave-ryzen/omega/imgview.html

// BOB Light Generator | CodePen
// https://codepen.io/NyteOwlDave/pen/ExJmPZQ

// Photorealism and Ray Tracing in C
// https://www.amazon.com/Photorealism-Ray-Tracing-Christopher-Watkins/dp/1558512470

// BOB for Windows (Ryzen)
// explorer "C:\home\dave\Apps\rt"
// file://C:/home/dave/Apps/rt

// NCS Graphics Depot (Tower)
// http://dave-tower/projects/gfx

// NCS Ray Tracers (Tower)
// http://dave-tower/projects/waypoints/raytracer

// MDN 2D Rendering Context
// https://developer.mozilla.org/en-US/docs/Web/API/CanvasRenderingContext2D/drawImage

// MDN Array Buffer
// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/ArrayBuffer

// Load Blob (see: snippets folder)
// https://chat.openai.com/c/18512bfa-28d7-4b07-ac7e-d079ef3dcca6

// Blob to Byte Array (see: snippets folder)
// https://chat.openai.com/c/55a65e50-9a90-42b9-80b7-7fd863bebee5

// Browser Page before unload and unload events
// https://chat.openai.com/c/5d94fc5b-f3f5-4667-a023-9d1903ba5090

// Platonic Solids
// http://dave-ryzen/examples/js/platonics.js

class BobException {
    constructor(object, member, details, category="none") {
        this.object = object;
        this.member = member;
        this.details = details;
        this.category = category;
    }
}

class BobToDoException extends BobException {
    constructor(object, member, details) {
        super(object, member, details, "TODO!");
    }
}

class BobTypeException extends BobException {
    constructor(object, member, argument, type) {
        super(object, member, 
            `Argument '${argument}' must be of type '${type}'`, 
            "Type Mismatch");
    }
}

const Bob = {
    precision: {
        MIN: 1, MAX: 6,
        value: 3
    },
    links: {
        home: "http://dave-tower/gfx/raytracer/bob"
    },
    visit: () => {
        window.open(Bob.links.home);
    },
    init: () => {
        try {
            console.group("BOB Ray Tracer Initialization");
        } catch(ex) {
            console.error(ex);
        } finally {
            console.groupEnd();
        }
    },
    shutdown: () => {
        try {
            console.group("BOB Ray Tracer Shutdown");
        } catch(ex) {
            console.error(ex);
        } finally {
            console.groupEnd();
        }
    },
    shutdownPending: (e, userConfirmCallback) => {
        try {
            console.group("BOB Ray Tracer Shutdown Requested");
            if ('function' == typeof userConfirmCallback) {
                Bob.userConfirmCallback = userConfirmCallback;
                return;
            } else {
                if ('function' == typeof Bob.userConfirmCallback) {
                    Bob.userConfirmCallback(e);
                }
            }
        } catch(ex) {
            console.error(ex);
        } finally {
            console.groupEnd();
        }
    },
    todo: (msg) => {
        console.warn("TODO", msg);
    },
};

window.addEventListener('load', e=>{
    Bob.init();
});

// Deprecated
// window.addEventListener('unload', e=>{
// });

window.addEventListener('beforeunload', e=>{
    Bob.shutdownPending( e );
});

class BobPoint {
    constructor(x, y=0) {
        this.acquire(x, y);
    }
    get primitive() {
        return {
            x: this.x,
            y: this.y
        };
    }
    get isNaN() {
        return isNaN(this.manhatten);
    }
    get isMajorX() {
        return this.w > this.h;
    }
    get isMajorY() {
        return this.w < this.h;        
    }
    get isDiagonal() {
        return this.x && (this.w == this.h);
    }
    get isOrigin() {
        return !(this.x || this.y);  
    }
    get x() {
        return this.__x__;
    }
    set x(value) {
        this.__x__ = this.fix(value);
    }
    get y() {
        return this.__y__;
    }
    set y(value) {
        this.__y__ = this.fix(value);        
    }
    get w() {
        return Math.abs(this.x);
    }
    get h() {
        return Math.abs(this.y);
    }
    get length() {
        return Math.hypot(this.x, this.y);
    }
    get lengthSquared() {
        return (this.x * this.x) + (this.y * this.y);
    }
    get manhatten() {
        return this.w + this.h;
    }
    get area() {
        return (this.w * this.h);        
    }
    get size() {
        return {
            w: this.w,
            h: this.h
        }
    }
    get rect() {
        return new BobRect(0, 0, this.x, this.y);
    }
    get perpendicular() {
        return new BobPoint(-this.y, this.x);
    }
    get angle() {
        return Math.atan2(this.y, this.x);
    }
    get aspect() {
        return this.w / this.h;
    }
    get slope() {
        return this.h / this.w;
    }
    get normal() {
        return this.vector.normal;
    }
    get vector() {
        return new BobVector(this.x, this.y, 0);
    }
    get quadrant() {
        if (this.isNaN ) return NaN;
        if (this.x > 0) {
            if (this.y > 0) {
                return 1; // Top-Right
            } else if (this.y < 0) {
                return 4; // Bottom-Right
            }    
        } else if (this.x < 0) {
            if (this.y > 0) {
                return 2;   // Top-Left
            } else if (this.y < 0) {
                return 3;   // Bottom-Left    
            }
        }
        return 0;
    }
    fix(n) {
        n = Number.parseInt(n);
        return isNaN(n) ? 0 : n;
    }
    clear() {
        this.__x__ = 
        this.__y__ = 0;
        return this;
    }
    clone() {
        return new BobPoint(this);
    }
    mimic(other) {
        if (other instanceof Object) {
            this.x = other.x;
            this.y = other.y;    
        } else {
            this.clear();
        }
        return this;
    }
    acquire(x, y) {
        if (x instanceof Object) {
            return this.mimic(x);
        } else {
            this.x = x;
            this.y = y;
            return this;
        }
    }
    flip() {
        return this.flipX().flipY();
    }
    flipX() {
        this.x = -this.x;
        return this;
    }
    flipY() {
        this.y = -this.y;
        return this;        
    }
    swapXY() {
        const temp = this.x;
        this.x = this.y;
        this.y = temp;
        return this;
    }
    add(other) {
        return new BobPoint(
            this.x + other.x,
            this.y + other.y
        );
    }
    sub(other) {
        return new BobPoint(
            this.x - other.x,
            this.y - other.y
        );
    }
    mul(other) {
        return new BobPoint(
            this.x * other.x,
            this.y * other.y
        );
    }
    dot(other) {
        return this.x * other.x + this.y * other.y;
    }
    dotPerp(other) {
        return this.x * other.y - this.y * other.x;
    }
    cross(other) {
        return new BobPoint(
              this.x * other.y,
            - this.y * other.x
        );
    }
    scale(k) {
        return new BobPoint(
            this.x * k,
            this.y * k
        );
    }
    lerp(t, other) {
        const v = this.vector.lerp(t, other.vector);
        return new BobPoint(v.x, v.y);
    }
    project(distance, direction) {
        const v = this.vector.project(distance, direction);
        return new BobPoint(v.x, v.y);
    }
    combine(ta, other, tb) {
        const v = this.vector.combine(ta, other.vector, tb);
        return new BobPoint(v.x, v.y);        
    }
    directionTo(other) {
        return other.sub(this).normal;
    }
    distanceTo(other) {
        return this.sub(other).length;
    }
    squaredDistanceTo(other) {
        return this.sub(other).lengthSquared;
    }
    getNearestIndex(points) {
        let minimum = Infinity;
        let index = -1;
        const limit = points.length;
        for (let i=0; i<limit; i++) {
            const dist = this.squaredDistanceTo(points[i]);
            if (dist < minimum) {
                minimum = dist;
                index = i;
            }
        }
        return index;
    }
    draw(ctx, color, radius=1, solid=true, border=true) {
        if (color) {
            const oldf = ctx.fillColor;
            const olds = ctx.strokeColor;
            ctx.strokeColor = ctx.fillColor = color;
            this.draw(ctx, radius);
            ctx.fillColor   = oldf;
            ctx.strokeColor = olds;    
        } else {
            // https://developer.mozilla.org/en-US/play
            ctx.beginPath();
            ctx.ellipse(this.x, this.y, radius, radius);
            if (solid ) ctx.fill();
            if (border) ctx.stroke();
        }
        return this;
    }
    drawSolid(ctx, color, radius=1) {
        return this.draw(ctx, color, radius, true, false);
    }
    drawBorder(ctx, color, radius=1) {
        return this.draw(ctx, color, radius, false, true);
    }
    toString() {
        return JSON.stringify(this.primitive);
    }
    parse(json) {
        return this.mimic(JSON.parse(json));
    }
}

class BobRect {
    constructor(left, top=0, right=0, bottom=0) {
        this.aquire(left, top, right, bottom)
    }
    get isNaN() {
        return isNaN(this.manhatten);
    }
    get isEmpty() {
        return this.area == 0;
    }
    get isNormal() {
        return (
            (this.right >= this.left) &&
            (this.bottom >= this.top)
        );
    }
    get isSquare() {
        return this.width == this.height;
    }
    get isWide() {
        return this.width > this.height;
    }
    get isTall() {
        return this.width < this.height;
    }
    set left(value) {
        const fix = this.fix;
        this.__left__ = fix(value);
    }
    set top(value) {
        const fix = this.fix;
        this.__top__ = fix(value);
    }
    set right(value) {
        const fix = this.fix;
        this.__right__ = fix(value);
    }
    set bottom(value) {
        const fix = this.fix;
        this.__bottom__ = fix(value);
    }
    get left() {
        return this.__left__;
    }
    get top() {
        return this.__top__;
    }
    get right() {
        return this.__right__;
    }
    get bottom() {
        return this.__bottom__;
    }
    get width() {
        return Math.abs(this.right - this.left);
    }
    get height() {
        return Math.abs(this.bottom - this.top);
    }
    get aspect() {
        return this.width / this.height;
    }
    get slope() {
        return this.height / this.width;
    }
    get size() {
        return {
            w: this.width,
            h: this.height
        };
    }
    get area() {
        return this.width * this.height;
    }
    get diagonal() {
        return Math.hypot(this.width, this.height);
    }
    get manhatten() {
        return this.width + this.height;
    }
    get perimeter() {
        return this.manhatten * 2;
    }
    get normal() {
        return new BobRect(this).normalize();
    }
    get leftTop() {
        return new BobPoint(this.left, this.top);
    }
    get leftBottom() {
        return new BobPoint(this.left, this.bottom);
    }
    get rightTop() {
        return new BobPoint(this.right, this.top);
    }
    get rightBottom() {
        return new BobPoint(this.right, this.bottom);
    }
    get center() {
        const x = (this.left + this.right) / 2;
        const y = (this.top + this.bottom) / 2;
        return new BobPoint(x, y);
    }
    get primitive() {
        return {
            left:   this.left,            
            top:    this.top,            
            right:  this.right,            
            bottom: this.bottom,            
        };
    }
    // Note: Y-axis is inverted!
    get quadrantIV() {
        const c = this.center;
        return new BobRect(c.x, this.top, this.right, c.y);
    }
    // Note: Y-axis is inverted!
    get quadrantIII() {
        const c = this.center;
        return new BobRect(this.left, this.top, c.x, c.y);
    }
    // Note: Y-axis is inverted!
    get quadrantII() {
        const c = this.center;
        return new BobRect(this.left, c.y, c.x, this.bottom);
    }
    // Note: Y-axis is inverted!
    get quadrantI() {
        const c = this.center;
        return new BobRect(c.x, c.y, this.right, this.bottom);
    }
    // Note: Y-axis is inverted!
    get quadrant() {
        if (this.x > 0) {
            return (this.y > 0) ? 1 : 4;
        } else if (this.x < 0) {
            return (this.y > 0) ? 2 : 3;
        }
        return 0;
    }
    fix(n) {
        n = Number.parseInt(n);
        return isNaN(n) ? 0 : n;
    }
    mimic(other) {
        if (other instanceof Object) {
            this.left   = other.left;
            this.right  = other.right;
            this.top    = other.top;
            this.bottom = other.bottom;    
        } else {
            this.clear();
        }
        return this;
    }
    acquire(left, top=0, right=0, bottom=0) {
        if (left instanceof Object) {
            return this.mimic(left);
        } 
        if (!this.isNaN(left)) {
            this.left   = left;
            this.right  = right;
            this.top    = top;
            this.bottom = bottom;
        } else {
            this.clear();
        }
        return this;
    }
    clear() {
        this.__left__   = 
        this.__right__  = 
        this.__top__    = 
        this.__bottom__ = 0;
        return this;
    }
    clone() {
        return new BobRect(this);
    }
    normalize() {
        const min = Math.min;
        const max = Math.max;
        // const mid = (a, b, c) => min(max(a, b), c);
        const x1 = this.left;
        const y1 = this.top;
        const x2 = this.right;
        const y2 = this.bottom;
        this.left   = min(x1, x2);
        this.top    = min(y1, y2);
        this.right  = max(x1, x2);
        this.bottom = max(y1, y2);
        return this;
    }
    flipX() {
        const temp = this.left;
        this.left = this.right;
        this.right = temp;
        return this;
    }
    flipY() {
        const temp = this.top;
        this.top = this.bottom;
        this.bottom = temp;
        return this;
    }
    flip() {
        return this.flipX().flipY();
    }
    swapXY() {
        let temp = this.left;
        this.left = this.top;
        this.top = temp;
        temp = this.right;
        this.right = this.bottom;
        this.bottom = temp;
        return this;
    }
    offset(dx, dy) {
        this.left   += dx;
        this.right  += dx;
        this.top    += dy;
        this.bottom += dy;
        return this;
    }
    inflate(dx, dy) {
        this.left   -= dx;
        this.right  += dx;
        this.top    -= dy;
        this.bottom += dy;
        return this;
    }
    intersect(other) {
        const max = Math.max;
        const min = Math.min;
        const a = this.normal;
        const b = other.normal;
        // Use quadrant or isNormal to detect non-intersection
        return new BobRect(
            max(a.left, b.left),
            max(a.top, b.top),
            min(a.right, b.right),    
            min(a.bottom, b.bottom)    
        );
    }
    union(other) {
        const max = Math.max;
        const min = Math.min;
        const a = this.normal;
        const b = other.normal;
        // Use isNaN or isEmpty to detect degenerate cases
        return new BobRect(
            min(a.left, b.left),
            min(a.top, b.top),
            max(a.right, b.right),    
            max(a.bottom, b.bottom)    
        );
    }
    equals(other) {
        return (
            (this.left   == other.left  ) &&
            (this.top    == other.top   ) &&
            (this.right  == other.right ) &&
            (this.bottom == other.bottom)
        );
    }
    contains(other) {
        return (
            (this.left   <= other.left  ) &&
            (this.top    <= other.top   ) &&
            (this.right  >= other.right ) &&
            (this.bottom >= other.bottom)
        );
    }
    containsPoint(pt) {
        return (
            (this.left   <= pt.x) &&
            (this.top    <= pt.y) &&
            (this.right  >= pt.x) &&
            (this.bottom >= pt.y)
        );
    }
    createElement(tagName) {
        const e = document.createElement(tagName);
        e.style.left     = this.left;
        e.style.top      = this.top;
        e.style.right    = this.right;
        e.style.bottom   = this.bottom;
        e.style.position = "absolute";
        return e;
    }
    createImageData() {
        return new ImageData(this.width, this.height);
    }
    createCanvas() {
        const e  = document.createElement("canvas");
        e.width  = this.width;
        e.height = this.height;
        return e;
    }
    createImage() {
        return new BobImage(this.width, this.height);
    }
    createRasterMap() {
        return new BobRasterMap(this.width, this.height);
    }
    draw(ctx, strokeColor, fillColor) {
        ctx.strokeColor = strokeColor || ctx.strokeColor;
        ctx.fillColor = fillColor || ctx.fillColor;
        ctx.beginPath();
        ctx.moveTo(this.left, this.top);
        ctx.lineTo(this.right, this.top);
        ctx.lineTo(this.right, this.bottom);
        ctx.lineTo(this.left, this.bottom);
        ctx.closePath();
        if (fillColor) ctx.fill();
        if (strokeColor) ctx.stroke();
        return this;
    }
    drawSolid(ctx, color) {
        return this.draw(ctx, null, color);
    }
    drawBorder(ctx, color) {
        return this.draw(ctx, color, null);
    }
    toString() {
        return JSON.stringify(this.primitive);
    }
    parse(json) {
        return this.mimic(JSON.parse(json));
    }
}

BobRect.fromCanvas = (o) => {
    return new BobRect(0, 0, o.width, o.height);
}
BobRect.fromElement = (o) => {
    o = o.style;
    return new BobRect(
        o.left, 
        o.top, 
        o.right, 
        o.bottom
    );
}
BobRect.fromImageData = (o) => {
    return new BobRect(0, 0, o.width, o.height);
}
BobRect.fromImage = (o) => {
    return new BobRect(0, 0, o.width, o.height);
}
BobRect.fromRasterMap = (o) => {
    return new BobRect(0, 0, o.width, o.height);
}

class BobVector {
    constructor(x, y=0, z=0) {
        this.acquire(x, y, z);
    }
    get primitive() {
        return {
            x: this.x,
            y: this.y,
            z: this.z
        };
    }
    get length() {
        return Math.hypot(this.x, this.y, this.z);
    }
    get lengthSquared() {
        return this.dot(this);
    }
    get manhatten() {
        const abs = Math.abs;
        return (
            abs(x) +
            abs(y) +
            abs(z)
        );
    }
    get normal() {
        return this.clone().normalize();
    }
    // Perspective projection into 2D space
    // on the 3D plane at Z=1
    // Points behind or at the eye (Z=0)
    // return a NULL vector 
    get flattened() {
        if (this.z > 0) {
            const k = 1 / this.z;
            return new BobVector(
                this.x * k,
                this.y * k,
                1
            );
        } else {
            return new BobVector();
        }
    }
    clear() {
        this.x =
        this.y =
        this.z = 0;
        return this;
    }
    clone() {
        return new BobVector(this);
    }
    mimic(other) {
        if (other instanceof Object) {
            this.x = other.x;
            this.y = other.y;
            this.z = other.z;    
        } else {
            this.clear();
        }
        return this;
    }
    acquire(x, y=0, z=0) {
        if (x instanceof Object) {
            return this.mimic(x);
        }
        if (!isNaN(x)) {
            this.x = x;
            this.y = y;
            this.z = z;    
        } else {
            this.clear();
        }
        return this;
    }
    add(other) {
        return new BobVector(
            this.x + other.x,
            this.y + other.y,
            this.z + other.z
        );
    }
    sub(other) {
        return new BobVector(
            this.x - other.x,
            this.y - other.y,
            this.z - other.z
        );
    }
    mul(other) {
        return new BobVector(
            this.x * other.x,
            this.y * other.y,
            this.z * other.z
        );
    }
    scale(k) {
        return new BobVector(
            this.x * k,
            this.y * k,
            this.z * k
        );
    }
    normalize() {
        const k = this.lengthSquared;
        if (k < 1e-08) {
            this.acquire(1);
        } else {
            const t = 1 / Math.sqrt(k);
            this.x *= t;
            this.y *= t;
            this.z *= t;
        }
        return this;
    }
    lerp(t, other) {
        const s = 1-t;
        return new BobVector(
            this.x * s + other.x * t,
            this.y * s + other.y * t,
            this.z * s + other.z * t
        );
    }
    project(t, other) {
        return new BobVector(
            this.x + other.x * t,
            this.y + other.y * t,
            this.z + other.z * t
        );
    }
    combine(ta, other, tb) {
        const s = 1-t;
        return new BobVector(
            this.x * ta + other.x * tb,
            this.y * ta + other.y * tb,
            this.z * ta + other.z * tb
        );
    }
    dot(other) {
        return (
            this.x * other.x + 
            this.y * other.y + 
            this.z * other.z
        );
    }
    cross(other) {
        return new BobVector(
            this.y * other.z - this.z * other.y,
            this.z * other.x - this.x * other.z,
            this.x * other.y - this.y * other.x
        );
    }
    gamma(n) {
        n = 1 / n;
        const g = x => Math.pow(x, n);
        return new BobVector(
            g(this.x),
            g(this.y),
            g(this.z)
        );
    }
    reflect(other) {
        const i = this.normal;      // I = Incident ray
        const n = other.normal;     // N = Surface normal
        const k = -i.dot(n) * 2;    // K = Scale
        return i.sub(n.scale(k));   // R = I - 2*(I dot N)*N
    }
    refract(other, ior) {
        // https://chat.openai.com/c/57bcbfdd-cfb5-43c5-9f37-6a5339e25892
        // https://en.wikipedia.org/wiki/List_of_refractive_indices
        ior = 1 / ior;                  // Simplification
        const i = this.normal;          // I = Incident ray
        const n = other.normal;         // N = Surface normal
        const klhs = i.dot(n);          // Klhs = Scale Left Hand Side
        const root = 1 - ior * ior * (1 - klhs * klhs);
        const krhs = Math.sqrt(root);   // Krhs = Scale Right Hand Side
        const lhs = i.sub(n.scale(k)).scale(ior);
        const rhs = n.scale(krhs);
        return lhs.sub(rhs);
    }
    directionTo(other) {
        return other.sub(this).normal;
    }
    distanceTo(other) {
        return other.sub(this).length;
    }
    squaredDistanceTo(other) {
        return other.sub(this).lengthSquared;
    }
    getNearestIndex(points) {
        let minimum = Infinity;
        let index = -1;
        const limit = points.length;
        for (let i=0; i<limit; i++) {
            const dist = this.squaredDistanceTo(points[i]);
            if (dist < minimum) {
                minimum = dist;
                index = i;
            }
        }
        return index;
    }
    toString() {
        return JSON.stringify(this.primitive);
    }
    parse(json) {
        return this.mimic(JSON.parse(json));
    }
}

class BobMatrix {
    constructor(other) {
        if (other instanceof Object) {
            this.clear();
        }
        this.mimic(other);
    }
    clear() {
        this.data = BobMatrix.create();
        return this;
    }
    clone() {
        return new BobMatrix(this);
    }
    mimic(other) {
        if (other instanceof Object) {
            BobMatrix.copy(other.data, this.data);
        } else {
            this.clear();
        }
        return this;
    }
    fill(n) {
        BobMatrix.fill(this.data, n);
        return this;
    }
    identity(n) {
        BobMatrix.identity(this.data);
        return this;
    }
    transpose() {
        const other = this.data;
        this.clear();
        BobMatrix.transpose(other, this.data);
        return this;
    }
    concat(other) {
        const result = new BobMatrix();
        BobMatrix.concat(this.data, other.data, result.data);
        return result;
    }
    apply(vin) {
        const vout = new BobVector();
        BobMatrix.concat(this.data, vin, vout);
        return vout;
    }
    toString() {
        return JSON.stringify(this.data);
    }
    parse(json) {
        const min = JSON.parse(json);
        BobMatrix.copy(min, this.data);
        return this;
    }
}

BobMatrix.create = () => {
    return [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]];
}
BobMatrix.fill = (mtx, n=0) => {
    for (let i=0; i<4; i++)
        for (let j=0; j<4; j++)
            mtx[i][j] = n;
}
BobMatrix.identity = (mtx) => {
    for (let i=0; i<4; i++)
        for (let j=0; j<4; j++)
            mtx[i][j] = (i==j) ? 1 : 0;
}
BobMatrix.copy = (min, mout) => {
    for (let i=0; i<4; i++)
        for (let j=0; j<4; j++)
            mout[i][j] = min[i][j];
}
BobMatrix.transpose = (min, mout) => {
    for (let i=0; i<4; i++)
        for (let j=0; j<4; j++)
            mout[i][j] = min[j][i];
}
BobMatrix.concat = (ma, mb, mout) => {
    BobMatrix.fill(mout);
    for (let i=0; i<4; i++)
        for (let j=0; j<4; j++)
            for (let k=0; k<4; kj++)
                mout[i][k] += ma[i][j] = mb[j][k];
}
BobMatrix.apply = (min, vin, vout) => {
    vout.x 
    = min[0][0] * vin.x
    + min[0][1] * vin.y
    + min[0][2] * vin.z
    + min[0][3];
    vout.y 
    = min[1][0] * vin.x
    + min[1][1] * vin.y
    + min[1][2] * vin.z
    + min[1][3];
    vout.z 
    = min[2][0] * vin.x
    + min[2][1] * vin.y
    + min[2][2] * vin.z
    + min[2][3];
}
BobMatrix.orthogonalize = (min, mout) => {
    Bob.todo("BobMatrix.orthogonalize");
}

class BobQuaternion {
    constructor(x, y=0, z=0, w=0) {
        this.acquire(x, y, z, w);
    }
    get primitive() {
        return {
            x: this.x,
            y: this.y,
            z: this.z,
            w: this.w
        };
    }
    get length() {
        return Math.hypot(
            this.x,
            this.y,
            this.z,
            this.w,
        );
    }
    get lengthSquared() {
        return this.dot(this);
    }
    get manhatten() {
        const abs = Math.abs;
        return (
            abs(x) +
            abs(y) +
            abs(z) +
            abs(w)
        );
    }
    get normal() {
        return this.clone().normalize();
    }
    // Perspective project into 3D space at the 4D volume W=1
    // Null and inverted volumes return a NULL quaternion
    get flattened() {
        if (this.w > 0) {
            const k = 1 / this.w;
            return new BobQuaternion(
                this.x * k,
                this.y * k,
                this.z * k,
                1
            );
        } else {
            return new BobQuaternion();
        }
    }
    clear() {
        this.x =
        this.y =
        this.z =
        this.w = 0;
        return this;
    }
    mimic(other) {
        if (other instanceof Object) {
            this.x = other.x;
            this.y = other.y;
            this.z = other.z;
            this.w = other.w;
        } else {
            this.clear();
        }
        return this;
    }
    acquire(x, y=0, z=0, w=0) {
        if (x instanceof Object) {
            return this.mimic(x);
        } 
        if (!isNaN(x)) {
            this.x = x;    
            this.y = y;    
            this.z = z;    
            this.w = w;    
        } else {
            this.clear();
        }
        return this;
    }
    clone() {
        return new BobQuaternion(this);
    }
    dot(other) {
        return (
            (this.x * other.x) +
            (this.y * other.y) +
            (this.z * other.z) +
            (this.w * other.w)
        );
    }
    add(other) {
        return new BobQuaternion(
            (this.x + other.x),
            (this.y + other.y),
            (this.z + other.z),
            (this.w + other.w)
        );
    }
    sub(other) {
        return new BobQuaternion(
            (this.x - other.x),
            (this.y - other.y),
            (this.z - other.z),
            (this.w - other.w)
        );
    }
    mul(other) {
        return new BobQuaternion(
            (this.x * other.x),
            (this.y * other.y),
            (this.z * other.z),
            (this.w * other.w)
        );
    }
    scale(k) {
        return new BobQuaternion(
            (this.x * k),
            (this.y * k),
            (this.z * k),
            (this.w * 
                k)
        );    
    }
    gamma(n) {
        n = 1 / n;
        const g = x => Math.pow(x, n);
        return new BobVector(
            g(this.x),
            g(this.y),
            g(this.z),
            g(this.w)
        );
    }
    lerp(t, other) {
        const s = 1-t;
        return new BobQuaternion(
            (this.x * s) + (other.x * t),
            (this.y * s) + (other.y * t),
            (this.z * s) + (other.z * t),
            (this.w * s) + (other.w * t)
        );
    }
    project(t, other) {
        return new BobQuaternion(
            this.x + (other.x * t),
            this.y + (other.y * t),
            this.z + (other.z * t),
            this.w + (other.w * t)
        );
    }
    combine(ta, other, tb) {
        return new BobQuaternion(
            (this.x * ta) + (other.x * tb),
            (this.y * ta) + (other.y * tb),
            (this.z * ta) + (other.z * tb),
            (this.w * ta) + (other.w * tb)
        );
    }
    normalize() {
        const k = this.lengthSquared;
        if (k < 1e-08) {
            this.acquire(1);
        } else {
            const t = 1 / Math.sqrt(k);
            this.x *= t;
            this.y *= t;
            this.z *= t;
            this.w *= t;
        }
        return this;
    }
    directionTo(other) {
        return other.sub(this).normal;
    }
    distanceTo(other) {
        return this.sub(other).length;
    }
    squaredDistanceTo(other) {
        return this.sub(other).lengthSquared;
    }
    getNearestIndex(points) {
        let minimum = Infinity;
        let index = -1;
        const limit = points.length;
        for (let i=0; i<limit; i++) {
            const dist = this.squaredDistanceTo(points[i]);
            if (dist < minimum) {
                minimum = dist;
                index = i;
            }
        }
        return index;
    }
    toString() {
        return JSON.stringify(this.primitive);
    }
    parse(json) {
        return this.mimic(JSON.parse(json));
    }
}

class BobColor {
    constructor(r, g=0, b=0) {
        this.acquire(r, g, b);
    }
    get primitive() {
        return {
            r: this.r,
            g: this.g,
            b: this.b
        };
    }
    get manhatten() {
        const abs = Math.abs;
        return abs(this.r) + abs(this.g) + abs(this.b);
    }
    get vector() {
        return new BobVector(
            this.r, 
            this.g, 
            this.b
        );
    }
    get isValid() {
        if (isNaN(this.r + this.g + this.b)) return false;
        if (this.r <    0) return false;
        if (this.g <    0) return false;
        if (this.b >    0) return false;
        if (this.r >  1.0) return false;
        if (this.g >  1.0) return false;
        if (this.b >  1.0) return false;
        return true;
    }
    clear() {
        this.r = 
        this.g = 
        this.b = 0;
        return this;
    }
    clone() {
        return new BobColor(this);
    }
    mimic(other) {
        if (other instanceof Object) {
            this.r = other.r;
            this.g = other.g;
            this.b = other.b;    
            return this.clamp();
        } else {
            this.clear();
            return this;
        }
    }
    acquire(r, g=0, b=0) {
        if (r instanceof Object) {
            return this.mimic(r);
        } 
        if (!isNaN(r)) {
            this.r = r;
            this.g = g;
            this.b = b;    
            return this.clamp();
        } else {
            return this.clear();
        }
    }
    clamp() {
        function fix(n) {
            if (isNaN(n)) return 0;
            return Math.min(Math.max(n, 0.0), 1.0);
        }
        this.r = fix(this.r); 
        this.g = fix(this.g); 
        this.b = fix(this.b);
        return this; 
    }
    toString() {
        return JSON.stringify(this.primitive);
    }
    parse(json) {
        return this.mimic(JSON.parse(json));
    }
}

BobColor.fromVector = (vec) => {
    return new BobColor(vec.x, vec.y, vec.z);
}

class BobPixel {
    constructor(r, g=0, b=0, a=255) {
        this.clear();
        this.acquire(r, g, b, a);
    }
    get primitive() {
        return {
            r: this.r,
            g: this.g,
            b: this.b,
            a: this.a
        };
    }
    get manhatten() {
        return this.data.reduce((a,b)=>a+b);
    }
    get r()     { return this.data[0]; }
    get g()     { return this.data[1]; }
    get b()     { return this.data[2]; }
    get a()     { return this.data[3]; }
    set r(n)    { this.data[0] = n; }
    set g(n)    { this.data[1] = n; }
    set b(n)    { this.data[2] = n; }
    set a(n)    { this.data[3] = n; }
    get bytes() {
        return new Uint8ClampedArray([
            this.r,
            this.g,
            this.b,
            this.a
        ]);
    }
    get integer() {
        return (
            (this.a << 24) |
            (this.r << 16) | 
            (this.g <<  8) | 
            this.b
        );
    }
    get color() {
        this.clamp();
        const frac = n => n / 255;
        return new BobColor( 
            frac(this.r),
            frac(this.g),
            frac(this.b)
        );
    }
    clear() {
        this.data = new Uint8ClampedArray([0,0,0,0]);
        return this;
    }
    clone() {
        return new BobPixel(this);
    }
    mimic(other) {
        if (Array.isArray(other)) {
            this.r = other[0];
            this.g = other[1];
            this.b = other[2];
            this.a = other[3];   
        } else if (other instanceof Object) {
            this.r = other.r;
            this.g = other.g;
            this.b = other.b;
            this.a = other.a;
        } else {
            this.clear();
        }
        return this;
    }
    acquire(r, g=0, b=0, a=255) {
        if (r instanceof Object) {
            return this.mimic(r);
        } else if (!isNaN(r)) {
            this.r = r;
            this.g = g;
            this.b = b;
            this.a = a;
        } else {
            this.clear();
        }
        return this;
    }
    toString() {
        return JSON.stringify(this.primitive);
    }
    parse(json) {
        return this.mimic(JSON.parse(json));
    }
}

BobPixel.fromColor = (color) => {
    const fix = n => {
        n = Math.round(n * 255);
        return isNaN(n) ? 0 : n;
    }
    color.clamp();
    return new BobPixel(
        fix(color.r),
        fix(color.g),
        fix(color.b)
    );
} 

BobPixel.fromInteger = (n) => {
    const a = (n >> 24) & 0xFF;
    const r = (n >> 16) & 0xFF;
    const g = (n >>  8) & 0xFF;
    const b = (n      ) & 0xFF;
    return new BobPixel(r, g, b, a);
}

BobPixel.fromBytes = (arr) => {
    return new BobPixel(
        arr[2],
        arr[1],
        arr[0],
        arr[3]
    );    
}

// Wrapper for ImageData
// https://hacks.mozilla.org/2011/12/faster-canvas-pixel-manipulation-with-typed-arrays/
class BobRasterMap {
    constructor(w, h) {
        if (w instanceof Object) {
            this.acquire(w);
        }  else {
            this.prepare(w || 1, h || 1);
        }
    }
    get width() {
        return this.image.width;
    }
    get height() {
        return this.image.height;
    }
    prepare(w, h) {
        this.image = new ImageData(w, h)
        return this;
    }
    clear() {
        return this.prepare(this.width, this.height);
    }
    clone() {
        return new BobRasterMap(this);
    }
    capture(image) {
        if (image instanceof ImageData) {
            const w = image.width;
            const h = image.height;
            const count = 4 * w * h;
            this.prepare(w, h);
            const src = image.data;
            const dst = this.image.data;
            for (let i=0; i<count; i++) {
                dst[i] = src[i];
            }
            return this;
        }
        throw new BobTypeException(
            "BobRasterMap",
            "capture", "image",
            "ImageData"
        );
    } 
    mimic(other) {
        if (other instanceof BobRasterMap) {
            return this.capture(other.image);
        }
        throw new BobTypeException(
            "BobRasterMap",
            "mimic", "other",
            "BobRasterMap"
        );
    }
    acquire(source) {
        if (source instanceof ImageData) {
            return this.capture(other);
        }
        if (source instanceof BobRasterMap) {
            return this.mimic(other);
        }
        this.image = BobRasterMap.acquire(source);
        return this;
    }
    getPixel(x, y) {
        const pixel = new BobPixel();
        return this.getPixelEx(x, y, pixel);
    }
    getPixelEx(x, y, pixel) {
        const img = this.image;
        const map = img.data;
        const index = 4*(y*img.width+x);
        pixel.b = map[index][0];
        pixel.g = map[index][1];
        pixel.r = map[index][2];
        pixel.a = map[index][3];
        return pixel;
    }
    putPixel(x, y, pixel) {
        const img = this.image;
        const map = img.data;
        const index = 4*(y*img.width+x);
        map[index][0] = pixel.b;
        map[index][1] = pixel.g;
        map[index][2] = pixel.r;
        map[index][3] = pixel.a;
        return this;
    }
    putPixelEx(x, y, r, g, b, a) {
        const img = this.image;
        const map = img.data;
        const index = 4*(y*img.width+x);
        map[index][0] = b;
        map[index][1] = g;
        map[index][2] = r;
        map[index][3] = a;
        return this;
    }
    toCanvas() {
        return this.toContext().canvas;
    }
    toContext() {
        const w = this.width;
        const h = this.height;
        const canvas = document.createElement('canvas');
        canvas.width = w;
        canvas.height = h;
        const ctx = canvas.getContext('2d');
        ctx.putImageData(this.image, 0, 0, w, h);
        return ctx;
    }
}

BobRasterMap.isImageData = e => {
    return e instanceof ImageData;
}

BobRasterMap.has = (imgData, x, y) => {
    const me = BobRasterMap;
    if (isNaN(x + y)) return false;
    if (me.isImageData(imgData)) {
        if (x < 0) return false;
        if (y < 0) return false;
        if (x >= imgData.width ) return false;
        if (y >= imgData.height) return false;
        return true;
    }
    return false;
}

BobRasterMap.getPixel = (imgData, x, y) => {
    const me = BobRasterMap;
    if (me.has(imgData, x, y)) {
        const index = 4*(y*imgData.width+x);
        const map = imgData.data;
        return new BobPixel(
            map[index][2],
            map[index][1],
            map[index][0],
            map[index][3]
        );
    }
    return new BobPixel();
}

BobRasterMap.putPixel = (imgData, x, y, pixel, alpha=0) => {
    const me = BobRasterMap;
    if (BobPixel.isPixel(pixel)) {
        if (me.has(imgData, x, y)) {
            const map = imgData.data;
            const index = 4*(y*imgData.width+x);
            map[index][0] = pixel.b;
            map[index][1] = pixel.g;
            map[index][2] = pixel.r;
            map[index][3] = pixel.a;
        }
    }
    return me;
}

BobRasterMap.acquire = (source) => {
    const me = BobRasterMap;
    if (source instanceof ImageData) {
        return me.fromImageData(source).image;
    }
    if (source instanceof BobImage) {
        return me.fromBobImage(source).image;
    }
    if (source instanceof BobRasterMap) {
        return me.fromBobRasterMap(source).image;
    }
    if (source instanceof HTMLCanvasElement) {
        return me.fromCanvas(source).image;
    }
    if (source instanceof CanvasRenderingContext2D) {
        return me.fromContext(source).image;
    }
    throw new TypeError("BobRasterMap.acquire");
}

BobRasterMap.fromImageData = (image) => {
    if (image instanceof ImageData) {
        return new BobRasterMap(image);
    }
    throw new BobTypeException(
        "BobRasterMap",
        "fromImageData", 
        "image",
        "ImageData"
    );
}

BobRasterMap.fromBobImage = (image) => {
    if (source instanceof BobImage) {
        image = new BobImage();
        return BobCanvasMap.fromCanvas(image.canvas);
    }
    throw new BobTypeException(
        "BobRasterMap",
        "fromBobImage", 
        "image",
        "BobImage"
    );
}

BobRasterMap.fromBobRasterMap = (map) => {
    if (source instanceof BobImage) {
        return new BobRasterMap(map);
    }
    throw new BobTypeException(
        "BobRasterMap",
        "fromBobRasterMap", 
        "map",
        "BobRasterMap"
    );
}

BobRasterMap.fromCanvas = (canvas) => {
    if (source instanceof HTMLCanvasElement) {
        const ctx = canvas.getContext('2d');
        return BobRasterMap.fromContext(ctx);
    }
    throw new BobTypeException(
        "BobRasterMap",
        "fromCanvas", 
        "canvas",
        "HTMLCanvasElement"
    );
}

BobRasterMap.fromContext = (context) => {
    if (source instanceof CanvasRenderingContext2D) {
        const map = new BobRasterMap();
        const canvas = context.canvas;
        const w = canvas.width;
        const h = canvas.height;
        map.image = context.getImageData(0, 0, w, h);
        return map;
    }
    throw new BobTypeException(
        "BobRasterMap",
        "fromContext", 
        "image",
        "CanvasRenderingContext2D"
    );
}

class BobImage {
    constructor() {
        this.clear()
    }
    clear() {
        this.canvas = null;
    }
    get ready() {
        return BobImage.isCanvas(this.canvas);
    }
    get width() {
        return this.ready ? this.canvas.width : 0;
    }
    get height() {
        return this.ready ? this.canvas.height : 0;
    }
    get size() {
        return {
            w: this.width,
            h: this.height
        }
    }
    render(canvas, rc) {
        if (this.ready && BobImage.isCanvas(canvas)) {
            const sw = this.width;
            const sh = this.height;
            const dw = rc.width;
            const dh = rc.height;
            const dx = rc.left;
            const dy = rc.top;
            const ctx = canvas.getContext2d();
            ctx.drawImage(
                this.canvas, 
                 0, 0, 
                sw, sh, 
                dx, dy, 
                dw, dh
            );
        }
        return this;
    }
    load(file) {    // File instance
        function read(file) {
            return new Promise((resolve, reject) => {
                const reader = new FileReader();    
                reader.onload  = () => { resolve(reader.result); };    
                reader.onerror = () => { reject(reader.error);   };
                reader.readAsArrayBuffer(file);
            });
        }
        function report(err) {
            console.error({
                message: "Error loading file",
                file,
                error: err,
                source: this
            });
        }
        function accept(arr) {
            this.decode(arr);
        }
        read(file).then(accept).catch(report);
    }
    // Arg is ArrayBuffer
    decode(arr) {
        const sr = new BobStreamReader(arr);
    }
}

BobImage.isCanvas = e => {
    return e instanceof HTMLCanvasElement;
}

// Encodes an image into a byte array in BOB image file format
// The image is in ImageData form
class BobImageEncoder {
    constructor(source) {
        if (source instanceof ImageData) {
            this.source = source;
        } else {
            throw new BobTypeException(
                "BobImageEncoder", "constructor", 
                "source", "ImageData"
            );
        }
    }
    encode(writer) {
        const w = this.source.width;
        const h = this.source.height;
        const pixels = this.source.data;
        writer.putw(w);
        writer.putw(h);
        writer.putdw(0);
        writer.putw(24);
        let index = 0, count = 1;
        let ro = pixels[index++];
        let go = pixels[index++];
        let bo = pixels[index++];
        ++index;
        const flush = (r, g, b) => {
            writer.putch(count);
            writer.putch(bo);
            writer.putch(go);
            writer.putch(ro);
            count = 1;
            ro = r;
            go = g;
            bo = b;
        }
        for (let y=0; y<h; y++) {
            for (let x=1; x<w; x++) {
                let r = pixels[index++];
                let g = pixels[index++];
                let b = pixels[index++];
                ++index;
                if ((r-ro)||(g-go)||(b-bo)||(count>253)) {
                    flush(r, g, b);
                } else {
                    ++count;
                }
            }
            flush(0, 0, 0);
        }
    }
}

// Decodes a byte array in BOB image file format into an image
// The image is in ImageData form
class BobImageDecoder {
    constructor(source) {
        if (source instanceof BobStreamReader) {
            this.source = source;
        } else {
            throw new BobTypeException(
                "BobImageDecoder", "constructor", 
                "source", "BobStreamReader"
            );
        }
    }
    decode() {
        // const reader = this.source;
        const reader = new BobStreamReader();
        reader.rewind();
    }
}

// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Guide/Typed_arrays
class BobStreamCore {
    constructor() {
        this.bigIndian = true;
        this.clear();
    }
    clear() {
        this.buffer = null;
        return this;
    }
    createBuffer(size) {
        this.buffer = BobStreamCore.createBuffer(size);
        return this.rewind();
    }
    get ready() {
        return BobStreamCore.isBuffer(this.buffer);
    }
    get length() {
        return this.ready ? this.buffer.length : 0;
    }
    get index() {
        return this["(index)"];
    }
    set index(n) {
        n = Math.min(n || 0, this.length);
        this["(index)"] = n;
    }
    has(index) {
        if (this.ready) {
            if (index < 0) return false;
            return (index < this.length);
        }
        return false;
    }
    peek(index) {
        if (this.ready) {
            index = isNaN(index) 
            ? this.index 
            : parseInt(index);
            return this.buffer[index];
        } else {
            return undefined;
        }
    }
}

BobStreamCore.isBuffer = (buffer) => {
    return (buffer instanceof Uint8Array);
}
BobStreamCore.createBuffer = (size) => {
    // const me = BobStreamCore;
    if (isNaN(size) || (size < 1)) {
        throw new Error(`Invalid buffer size: ${size}`);
    }
    return new Uint8Array(size);
}


// Converts various object types into
// buffers of type Uint8Array (unsigned byte arrays)
BobStreamCore.acquireBuffer = (source) => {
    function todo(details) {
        throw new BobToDoException(
            'BobStreamCore',
            'acquireBuffer',
            details
        );
    }
    async function arr(s) {
        return await s.arrayBuffer();
    }
    function bytes(s) {
        const count = s.length;
        const result = new Uint8Array(length);
        for (let i=0; i<count; i++) {
            result[i] = s[i];
        }
        return result;
    }
    if (source instanceof Blob) {
        return bytes(arr(source));
    } else if (source instanceof ArrayBuffer) {
        return bytes(source);
    } else if (Array.isArray(source)) {
        return bytes(source);
    } else if (source instanceof ImageData) {
        return bytes(source.data);
    }
    throw new TypeError("BobStreamCore.acquireBuffer");
}
BobStreamCore.growBuffer = (buffer, addBytes) => {
    const me = BobStreamCore;
    const assert = (test, message) => {
        console.assert(test, {
            message, buffer, addBytes,
            source: "BobStreamCore.growBuffer"
        });    
    }
    const ok = Array.isArray;
    assert(addBytes >= 0, "Added byte count can't be < 0");
    if (!addBytes) return buffer; // Nothing to do?
    assert(ok(buffer), "Invalid buffer data type");
    let bufLen = buffer.length;
    let totalLen = bufLen + addBytes;
    const create = me.createBuffer;
    const pending = create(totalLen);
    for (let i=0; i<bufLen; i++) {
        pending[i] = buffer[i];
    }
    return pending;
}
BobStreamCore.resizeBuffer = (buffer, totalBytes) => {
    const me = BobStreamCore;
    const assert = (test, message) => {
        console.assert(test, {
            message, buffer, totalBytes,
            source: "BobStreamCore.resizeBuffer"
        });    
    }
    const ok = Array.isArray;
    assert(totalBytes > 0, "Total byte count must be > 0");
    assert(ok(buffer), "Invalid buffer data type");
    let bufLen = buffer.length;
    if (bufLen == totalBytes) return buffer;    // Nothing to do?
    if (bufLen > totalBytes) {                  // Need to grow?
        return me.growBuffer(buffer, totalBytes - bufLen); 
    }
    const create = me.createBuffer;
    const pending = create(totalBytes);
    for (let i=0; i<totalBytes; i++) {
        pending[i] = buffer[i];
    }
    return pending;
}
BobStreamCore.appendBuffer = (dst, src) => {
    const me = BobStreamCore;
    const assert = (test, message) => {
        console.assert(test, {
            message, 
            src, dst,
            source: "BobStreamCore.appendBuffer"
        });
    }
    const ok = Array.isArray;
    assert(ok(dst), "Invalid dst buffer data type");
    assert(ok(src), "Invalid src buffer data type");
    const pending = me.growBuffer(dst, src.length);
    for (let o=dst.length, i=0; i < src.length; i++) {
        pending[o+i] = src[i];
    }
    return pending;
}

// Reads a stream from a buffer
class BobStreamReader extends BobStreamCore {
    constructor(source) {
        super();
        this.acquireBuffer(source);
    }
    get bof() {
        if (super.ready) {
            return this.index < 1;
        } else {
            return true;
        }
    }
    get eof() {
        return this.remaining ? false : true;
    }
    get remaining() {
        if (super.ready) {
            const count = this.length;
            const index = this.index;
            if (index < 0) return 0;
            if (index >= count) return 0;
            return count - index;
        } else {
            return 0;
        }
    }
    clear() {
        super.clear();
        this.rewind();
        return this;
    }
    acquireBuffer(source) {
        if (source) {
            this.buffer = BobStreamCore.acquireBuffer(source);
            this.rewind();    
        } else {
            this.clear();
        }
        return this;
    }
    rewind() {
        this.index = 0;
        return this;
    }
    prev() {
        if (this.bof) return false;
        this.index = this.index - 1;
        return true;
    }
    next() {
        if (this.eof) return false;
        this.index = this.index + 1;
        return true;
    }
    getch() {
        const ch = this.peek();
        this.next();
        return ch;
    }
    ungetch(byte) {
        if (this.prev()) {
            this.buffer[this.index] = byte;
            return true;
        } else {
            return false;
        }
    }
    unget(bytes) {
        let count = 0;
        while (this.prev()) {
            if (!bytes.length) break;
            this.buffer[this.index] = bytes.shift();
            ++count;
        }
        return count;
    }
    getw() {
        const a = this.getch();
        const b = this.getch();
        if (this.bigEndian) {
            return (a << 8) | b;
        } else {
            return (b << 8) | b;
        }
    }
    getdw() {
        const a = this.getch();
        const b = this.getch();
        const c = this.getch();
        const d = this.getch();
        if (this.bigEndian) {
            return (a << 24) 
                 | (b << 16)
                 | (c <<  8)
                 | (d);
        } else {
            return (d << 24) 
                 | (c << 16)
                 | (b <<  8)
                 | (a);
        }
    }
    read(count) {
        count = Math.min(count || 0, this.remaining);
        if (count < 1) return 0;
        const bytes = new Uint8Array(count);
        for (let i=0; i<count; i++) {
            bytes[i] = this.buffer[this.index++];
        }
        return bytes;
    }
}

// Writes a stream to a buffer
class BobStreamWriter extends BobStreamCore {
    constructor(blocks=1) {
        super();
        this.createCache(blocks);
    }
    createCache(blocks=1) {
        blocks = parseInt(Math.max(1, blocks || 1));
        const bytes = cacheSize = blocks * BobStreamWriter.blockSize;
        this.cache = new Uint8Array(bytes);
        return this;
    }
    flush() {
        if (this.lengthCache < 1) return this;
        const append = BobStreamCore.appendBuffer;
        this.buffer = append(this.buffer, this.cache);
        return this.createCache(this.cacheBlocks);
    }
    flushIfCacheFull() {
        if (this.cacheRemaining < 1) return this.flush(); 
        return this;        
    }
    get cache() {
        return this["(cache)"];
    }
    set cache(buffer) {
        if (buffer instanceof Uint8Array) {
            this["(cache)"] = buffer;
        } else {
            throw new BobTypeException(
                "BobStreamWriter", "cache",
                "buffer", "Uint8Array"
            );
        }
    }
    get cacheIndex() {
        return this["(cacheIndex)"]
    }
    set cacheIndex(n) {
        n = Math.min(n || 0, this.cacheLength);
        this["(cacheIndex)"] = n;
    }
    get cacheLength() {
        return this.cache.length;
    }
    get cacheBlocks() {
        const bytes = this.cacheLength;
        const k = BobStreamWriter.blockSize;
        return (bytes / k) + ((bytes % k) ? 1 : 0);
    }
    get cacheRemaining() {
        return this.cacheLength - this.cacheIndex;
    }
    get lengthTotal() {
        return this.length + this.cacheLength;
    }
    putch(byte) {
        this.flushIfCacheFull();
        const cache = this.cache;
        const index = this.cacheIndex;
        cache[index++] = byte;
        this.cacheIndex = index;
        return this;
    }
    putw(word) {
        let a, b;
        if (this.bigEndian) {
            a = 0xFF & (word >> 8);
            b = 0xFF & word;
        } else {
            b = 0xFF & (word >> 8);
            a = 0xFF & word;
        }
        this.putch(a);
        this.putch(b);
    }
    putdw(dword) {
        let a, b, c, d;
        if (this.bigEndian) {
            a = 0xFF & (dword >> 24);
            b = 0xFF & (dword >> 16);
            c = 0xFF & (dword >>  8);
            d = 0xFF & dword;
        } else {
            d = 0xFF & (dword >> 24);
            c = 0xFF & (dword >> 16);
            b = 0xFF & (dword >>  8);
            a = 0xFF & dword;
        }
        this.putch(a);
        this.putch(b);
        this.putch(c);
        this.putch(d);
        return this;
    }
    write(bytes) {
        const count = bytes.length;
        for (let i=0; i<count; i++) {
            this.putch(bytes[i]);
        }
        return this;
    }
}

// Write cache block size
BobStreamWriter.blockSize = 512;

// === FUTURE ENHANCEMENTS 

// Ray
class BobRay {
    constructor() {
        Bob.todo("BobRay");
    }
}

// Plane
class BobPlane {
    constructor() {
        Bob.todo("BobPlane");
    }
}

// Light Source
class BobLight {
    constructor() {
        Bob.todo("BobLight");
    }
}

// Caustic Effects
class BobCaustics {
    constructor() {
        Bob.todo("BobCaustics");
    }
}

// Pixmap Texture
class BobTexMap {
    constructor() {
        Bob.todo("BobTexMap");
    }
}

// Procedural Texture
class BobTexProc {
    constructor() {
        Bob.todo("BobTexProc");
    }
}

// Either Type Texture
class BobTexture {
    constructor() {
        Bob.todo("BobTexture");
    }
}

// Collection of Textures
class BobTextureStack {
    constructor() {
        Bob.todo("BobTextureStack");
    }
}

// Height, Width, and Depth
class BobVolume {
    constructor() {
        Bob.todo("BobVolume");
    }
}

// Height and Width
class BobArea {
    constructor() {
        Bob.todo("BobArea");
    }
}

// Cuboid
class BobVoxel {
    constructor() {
        Bob.todo("BobVoxel");
    }
}

// View Frustum
class BobFrustrum {
    constructor() {
        Bob.todo("BobFrustum");
    }
}

// Tracer/Rasterizer
class BobTracer {
    constructor() {
        Bob.todo("BobTracer");
    }
}

// Studio Settings
class BobStudio {
    constructor() {
        Bob.todo("BobStudio");
    }
}

// Camera
class BobCamera {
    constructor() {
        Bob.todo("BobCamera");
    }
}

// Camera Lens
class BobCameraLens {
    constructor() {
        Bob.todo("BobCameraLens");
    }
}

// Viewpoint
class BobViewpoint {
    constructor() {
        Bob.todo("BobViewpoint");
    }
}

// Generic Shape
class BobShape {
    constructor() {
        Bob.todo("BobShape");
    }
}

// Polyhedron Shape
class BobPolyhedron {
    constructor() {
        Bob.todo("BobPolyhedron");
    }
}

// Ring Shape
class BobRing {
    constructor() {
        Bob.todo("BobRing");
    }
}

// Sphere Shape
class BobSphere {
    constructor() {
        Bob.todo("BobSphere");
    }
}

// Patch Shape
class BobPatch {
    constructor() {
        Bob.todo("BobPatch");
    }
}

// Polygon Shape
class BobPolygon {
    constructor() {
        Bob.todo("BobPolygon");
    }
}

// Clipping Volume
class BobClip {
    constructor() {
        Bob.todo("BobClip");
    }
}

// Perlin Noise
class BobNoise {
    constructor() {
        Bob.todo("BobNoise");
    }
}

// Transform
class BobTransform {
    constructor(vt, vs, vr) {
        this.vt = vt.clone();
        this.vs = vs.clone();
        this.vr = vr.clone();
    }
    createTSR() {
        Bob.todo("BobTransform.createTSR");
    }
    createTranslation() {
        Bob.todo("BobTransform.createTranslation");
    }
    createScaling() {
        Bob.todo("BobTransform.createScaling");
    }
    createRotation() {
        Bob.todo("BobTransform.createRotation");
    }
    createView() {
        Bob.todo("BobTransform.createView");
    }
}

// Transform Stack
class BobTransformStack {
    constructor() {
        Bob.todo("BobTransformStack");
    }
}

// Bounding Slabs Clipping
class BobBoundingSlabs {
    constructor() {
        Bob.todo("BobBoundingSlabs");
    }
}

// Color Palette
class BobPalette {
    constructor() {
        Bob.todo("BobPalette");
    }
}

// Color Palette Entry
class BobPaletteEntry {
    constructor(key, r, g, b) {
        this.key = key;
        this.r = r;
        this.g = g;
        this.b = b;
    }
    get html() {
        const hex = n => n.toString(16)[0].toUpperCase();
        this.clamp();
        if (this.isDim) {
            const r = hex(r);
            const g = hex(g);
            const b = hex(b);
            return `#${r}${g}${b}`;
        } else {
            let n = this.integer;
            let s = n.toString(16).toUpperCase();
            while (s.length < 6) {
                s = '0' + s;
            }
            return '#' + s;
        }
    }
    set html(text) {
        // TODO: Should we mess with abbreviated format here?
        // e.g. #FED => #F0E0D0
        if (text.startsWith('#')) {
            text = '0x' + text.substring(1);
        } else if (!(text.toLowerCase().startsWith('0x'))) {
            text = '0x' + text;
        }
        this.integer = Number.parseInt(text);
    }
    get css() {
        this.clamp();
        return `rgb(${this.r}, ${this.g}, ${this.b})`;
    }
    set css(text) {
        if (text.startsWith('rgb(')) {
            if (text.endsWith(')')) {
                const rgb = text.substring(4, text.length - 1)
                .split(',')
                .map(s=>s.trim());
                this.acquire(rgb[0], rgb[1], rgb[2]);
                return;
            }
        }
        this.clear();
    }
    get pixel() {
        return new BobPixel(
            this.r,
            this.g,
            this.b,
            0
        );
    }
    set pixel(other) {
        this.r = pixel.r;
        this.g = pixel.g;
        this.b = pixel.b;
    }
}

// Background
class BobBackground {
    constructor() {
        Bob.todo("BobBackground");
    }
}

