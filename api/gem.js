
// GEM was originally implement in the Playground
// http://dave-tower/projects/playground/
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// https://en.wikipedia.org/wiki/List_of_mathematical_shapes

function throw_not_implemented() {
    throw new Error("This feature is not yet implemented.");
}

function throw_out_of_range() {
    throw new Error("Argument out of range.");
}

function throw_matrix_size_mismatch() {
    throw new Error("Matrix size mismatch.");
}

function throw_vector_size_mismatch() {
    throw new Error("Vector size mismatch.");
}

// Another way to use typeof (prettier if you ask me)
function is_type(arg, type) {
    return (typeof arg === type);
}

// Whether the argument is iterable (for-in)
function is_iterable(arg) {
    if (arg) {
        return is_type(arg[Symbol.iterator], 'function');
    } else {
        return false;
    }
}

// Whether one of the argument's properties 
// is enumerable (for-of)
function is_enumerable_prop(arg, prop) {
    if (arg && is_type(arg, "object")) {
        return arg.propertyIsEnumerable(prop);
    } else {
        return false;
    }
}

// Logarithm to the base of n
exports.logn = function logn(n, x) {
    return Math.log(x) / Math.log(n);
}

// Nth root of x
exports.root = function root(x, n) {
    return Math.pow(x, 1/n);
}

// If n<0 throws an exception
exports.factorial = function factorial(n) {
    n = parseInt(n) || 0;
    if (n < 0) throw_out_of_range();
    if (n < 2) return 1;
    return n * factorial(n - 1);
}

// If n<0 throws an exception
exports.fibonacci = function fibonacci(n) {
    n = parseInt(n) || 0;
    if (n < 0) throw_out_of_range();
    if (n < 2) {
      return n;
    } else {
      return fibonacci(n - 1) + fibonacci(n - 2);
    }
}

exports.exp = Math.exp;
exports.expm1 = Math.expm1;
exports.log = Math.log;
exports.log10 = Math.log10;
exports.log1p = Math.log1p;
exports.log2 = Math.log2;
exports.pow = Math.pow;
exports.sqrt = Math.sqrt;
exports.cbrt = Math.cbrt;

// Cubes x
exports.cube = function cube(x) {
    return x*x*x;
}

// Squares x
exports.square = function square(x) {
    return x*x;
}

exports.sin = Math.sin;
exports.cos = Math.cos;
exports.tan = Math.tan;

exports.asin = Math.asin;
exports.acos = Math.acos;
exports.atan = Math.atan;
exports.atan2 = Math.atan2;

exports.abs = Math.abs;
exports.min = Math.min;
exports.max = Math.max;
exports.ceil = Math.ceil;
exports.floor = Math.floor;
exports.round = Math.round;
exports.trunc = Math.trunc;

exports.sinh = Math.sinh;
exports.cosh = Math.cosh;
exports.tanh = Math.tanh;

exports.asinh = Math.asinh;
exports.acosh = Math.acosh;
exports.atanh = Math.atanh;

// Handy constants
exports.pi = Math.PI;
exports.e = Math.E
exports.sqrt2 = Math.SQRT2;
exports.sqrt1_2 = Math.SQRT1_2;
exports.tiny = 1e-8;
exports.huge = 1e+8;

exports.clz32 = Math.clz32;
exports.fround = Math.fround;
exports.hypot = Math.hypot;
exports.imul = Math.imul;
exports.rnd = Math.random;
exports.sgn = Math.sign;

// Similar to sgn() but returns 0 when x is zero
// or NaN
exports.sgnz = function sgnz(x) {
    return x ? sgn(x) : 0;
}

// Performs gamma correction on a single color channel
exports.gamma = function gamma(n,lo,hi,gam) {
    const range = hi - lo;
    const value =  n - lo;
    return lo + Math.pow(value/range, 1/gam) * range;
}

// A is the origin
// B is the terminus
// t is the interpolant
// Generally, t should be in the range [0.0 ... 1.0]
// When t is outside this range, extrapolation results
exports.lerp = function lerp(a, b, t) {
    return a * (t-1) + b * t;
}

// Ray projection in 1D
// A is th origin
// B is the direction normal [-1.0 ... +1.0]
// t is the distance
exports.project = function project(A,B,t) {
    return A + t*B;
}

// Vector combination in 1D
// A and B are the direction normals [-1.0 ... +1.0]
// a and b are the corresponding distances
exports.combine = function combine(A,a,B,b) {
    return a*A + b*B;
}

// Convert degrees to radians
exports.deg2rad = function deg2rad(degrees) {
    return degrees*Math.PI/180;
}

// Convert radians to degrees
exports.rad2deg = function rad2deg(radians) {
    return radians*180/Math.PI;
}

// Find middle value
exports.mid = function clamp(n, lo, hi) {
    return Math.min(Math.max(n,lo),hi);
}
  
// Modulus with guarantee of positiveness
exports.wrap = function wrap(n, limit) {
    n %= limit;
    // The +0 prevents -0 (which javascript supports)
    return ((n < 0) ? (n + limit) : (n)) + 0;
}

// Convert number to hexadecimal string
exports.hex = function hex(n) { return '0x' + parseInt(n).toString(16); }

// Convert number to decimal integer string
exports.dec = function dec(n) { return parseInt(n).toString(); }

// Convert number to decimal floating point string
exports.flt = function flt(n) {return parseFloat(n).toString(); }

// Convert number to binary string
exports.bin = function bin(n) { return '0b' + parseInt(n).toString(2); }

// Convert number to octal string
exports.oct = function oct(n) { return '0o' + parseInt(n).toString(8); }

exports.Vec2 = class Vec2 {
    constructor(x=0, y=0) {
        this.x = x;
        this.y = y;
    }
    get zero() {
        return new Object(Vec2.zero());
    }
    get x_axis() {
        return new Object(Vec2.x_axis());
    }
    get y_axis() {
        return new Object(Vec2.y_axis());
    }
    get vector() {
        new Vector().acquire([this.x, this.y]);
    }
    get polar() {
        return new VecPolar().acquire(this);
    }
    get complex() {
        return new Complex().acquire(this);
    }
    get clone() {
        return new Vec2().acquire(this);
    }
    get mod() {
        return Math.hypot(this.x, this.y);
    }
    get arg() {
        return Math.atan2(this.y, this.x);        
    }
    get normal() {
        let k = this.mod;
        if (isNaN(k) || (k < exports.tiny)) {
            return new Vec2(1, 0);
        } else {
            k = 1 / k;
            return new Vec2(this.x*k, this.y*k);
        }
    }
    set normal(v) {
        v = v.normal;
        this.x = v.x;
        this.y = v.y; 
    }
    acquire(v) {
        if (!Array.isArray(v)) {
            if (v instanceof Vector) v = v.vector.v;
            else v = [];
        }
        this.x = v[0] || 0;
        this.y = v[1] || 0;
        return this;
    }
    // Perspective projection into 1-space (scalar)
    flatten(origin=0, aspect=1) {
        const k = 1 / this.y;
        return this.x * k * aspect + origin;
    }
}

exports.Vec2.zero   = () => [0, 0];
exports.Vec2.x_axis = () => [1, 0];
exports.Vec2.y_axis = () => [0, 1];

exports.Vec3 = class Vec3 {
    constructor(x=0, y=0, z=0) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
    get zero() {
        return new Object(Vec3.zero());
    }
    get x_axis() {
        return new Object(Vec3.x_axis());
    }
    get y_axis() {
        return new Object(Vec3.y_axis());
    }
    get z_axis() {
        return new Object(Vec3.y_axis());
    }
    get vector() {
        return new Vector().acquire([this.x, this.y, this.z]);
    }
    get vec2() {
        return new Vec2().acquire(this);
    }
    get veccyl() {
        return new VecCyl().acquire(this);
    }
    get vecsph() {
        return new VecSph().acquire(this);
    }
    get clone() {
        return new Vec3().acquire(this);
    }
    get mod() {
        return Math.hypot(this.x, this.y, this.z);
    }
    get arg() {
        const n = this.normal;
        return new DirectionAngles().cosines =
            [
                n.x,
                n.y,
                n.z
            ];
    }
    get normal() {
        let k = this.mod;
        if (isNaN(k) || (k < exports.tiny)) {
            return new Vec3(1, 0, 0);
        } else {
            k = 1 / k;
            return new Vec3(this.x*k, this.y*k, this.z*k);
        }
    }
    set normal(v) {
        v = v.normal;
        this.x = v.x;
        this.y = v.y; 
        this.z = v.z; 
    }
    acquire(v) {
        if (!Array.isArray(v)) {
            if (v instanceof Vector) v = v.vector.v;
            else v = [];
        }
        this.x = v[0] || 0;
        this.y = v[1] || 0;
        this.z = v[2] || 0;
        return this;
    }
    // Perspective projection into 2-space
    flatten(origin={x:0, y:0}, aspect=1) {
        const k = 1 / this.z;
        const x = this.x * k * aspect + origin.x;
        const y = this.y * k + origin.y;
        return new Vec2(x, y);
    }
}

exports.Vec3.zero   = () => [0, 0, 0];
exports.Vec3.x_axis = () => [1, 0, 0];
exports.Vec3.y_axis = () => [0, 1, 0];
exports.Vec3.z_axis = () => [0, 0, 1];

exports.Vec4 = class Vec4 {
    constructor(x=0, y=0, z=0, w=0) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }
    get zero() {
        return new Object(Vec4.zero());
    }
    get x_axis() {
        return new Object(Vec4.x_axis());
    }
    get y_axis() {
        return new Object(Vec4.y_axis());
    }
    get z_axis() {
        return new Object(Vec4.y_axis());
    }
    get w_axis() {
        return new Object(Vec4.w_axis());
    }
    get vector() {
        return new Vector().acquire([
            this.x,
            this.y,
            this.z,
            this.w
        ]);
    }
    get vec3() {
        return new Vec3().acquire(this);
    }
    get clone() {
        return new Vec4().acquire(this);
    }
    get mod() {
        return Math.hypot(this.x, this.y, this.z, this.w);
    }
    get arg() {
        const n = this.normal;
        return new DirectionAngles().cosines =
            [
                n.x,
                n.y,
                n.z,
                n.w
            ];
    }
    get normal() {
        let k = this.modulus;
        if (isNaN(k) || (k < exports.tiny)) {
            return new Vec4(1, 0, 0, 0);
        } else {
            k = 1 / k;
            return new Vec4(this.x*k, this.y*k, this.z*k, this.w*k);
        }
    }
    set normal(v) {
        v = v.normal;
        this.x = v.x;
        this.y = v.y; 
        this.z = v.z; 
        this.w = v.w; 
    }
    acquire(v) {
        if (!Array.isArray(v)) {
            if (v instanceof Vector) v = v.vector.v;
            else v = [];
        }
        this.x = v[0] || 0;
        this.y = v[1] || 0;
        this.z = v[2] || 0;
        this.w = v[3] || 0;
        return this;
    }
    // Perspective projection into 3-space
    flatten(origin={x:0, y:0, z:0, w:0}, aspect={x:1, y:1, z:1, w:1}) {
        const k = 1 / this.w;
        const x = this.x * k * aspect.x + origin.x;
        const y = this.y * k * aspect.y + origin.y;
        const z = this.z * k * aspect.z + origin.z;
        const w = this.w * k * aspect.w + origin.w;
        return new Vec3(x, y, z, w);
    }
}

exports.Vec4.zero   = () => [0, 0, 0, 0];
exports.Vec4.x_axis = () => [1, 0, 0, 0];
exports.Vec4.y_axis = () => [0, 1, 0, 0];
exports.Vec4.z_axis = () => [0, 0, 1, 0];
exports.Vec4.w_axis = () => [0, 0, 0, 1];

// https://en.wikipedia.org/wiki/Cylindrical_coordinate_system
exports.VecCyl = class VecCyl {
    constructor(rho=0, theta=0, height=0) {
        this.rho = rho;
        this.theta  = theta;
        this.height = height;
    }
    get vec3() {
        throw_not_implemented();
    }
    get vecsph() {
        return new VecSph().aquire(this.vec3);
    }
    get mod() {
        return this.vec3.mod;
    }
    get arg() {
        return this.vec3.arg;
    }
    get clone() {
        return new VecCyl().acquire(this.vec3);
    }
    get normal() {
        return new VecCyl.acquire(this.vec3.normal);
    }
    set normal(v) {
        return this.acquire(new Vec3.acquire(v).normal);
    }
    acquire(v) {
        /*
        v = new Vec3().acquire(v);
        const rho = 0;
        const theta = 0;
        const height = 0;
        this.rho = rho;
        this.theta = theta;
        this.height = height;
        return this;
        */
       throw_not_implemented();
    }
}

// https://en.wikipedia.org/wiki/Spherical_coordinate_system
exports.VecSph = class VecSph {
    constructor(rho=0, theta=0, phi=0) {
        this.rho = rho;
        this.theta  = theta;
        this.phi = phi;
    }
    get vec3() {
        throw_not_implemented();
    }
    get veccyl() {
        return new VecCyl().aquire(this.vec3);
    }
    get mod() {
        return this.rho;
    }
    get arg() {
        return this.vec3.arg;
    }
    get clone() {
        return new VecSph().acquire(this.vec3);
    }
    get normal() {
        return new VecSph.acquire(this.vec3.normal);
    }
    set normal(v) {
        return this.acquire(new Vec3.acquire(v).normal);
    }
    acquire(v) {
        throw_not_implemented();
        /*
        v = new Vec3().acquire(v);
        const rho = 0;
        const theta = 0;
        const phi = 0;
        this.rho = rho;
        this.theta = theta;
        this.phi = phi;
        return this;
        */
    }
}

// https://en.wikipedia.org/wiki/Polar_coordinate_system
exports.VecPolar = class VecPolar {
    constructor(mod=0, arg=0) {
        this.mod = mod;
        this.arg = arg;
    }
    get vec2() {
        const x = this.mod * Math.cos(this.arg);
        const y = this.mod * Math.sin(this.arg);
        return new Vec2(x, y);
    }
    get complex() {
        return new Complex().acquire(this.vec2);
    }
    // NOTE: mod and arg are properties
    get clone() {
        return new VecPolar(this.mod, this.arg);
    }
    get normal() {
        return new VecPolar().acquire(this.vec2.normal);
    }
    set normal(v) {
        return this.acquire(new Vec2.acquire(v).normal);
    }
    acquire(v) {
        v = new Vec2().acquire(v);
        this.mod = v.mod;
        this.arg = v.arg;
        return this;
    }
}

// https://en.wikipedia.org/wiki/Vector_(mathematics_and_physics)
exports.Vector = class Vector {
    constructor(other=null) {
        if (other) {
            this.acquire(other);
        } else {
            this.length = 0;
        }
    }
    get length() {
        return this.v.length;
    }
    set length(n) {
        this.v = new Array(n).fill(0);
    }
    get clone() {
        return new Vector(this);
    }
    get reciprocals() {
        const vout = new Vector();
        const count = vout.length = this.length;
        for(let i=0; i<count; i++) {
            vout.v[i] = 1 / this.v[i];
        }
        return vout;        
    }
    get sum() {
        return this.v.reduce((total, num)=>total+num);
    }
    get sum_of_squares() {
        return this.v.reduce((total, num)=>total+num*num);        
    }
    get modulus() {
        return Math.sqrt(this.sum_of_squares);
    }
    get mean() {
        if (this.length) {
            return this.sum_of_squares / this.length;
        }
        return 0;

    }
    get median() {
        const index = Math.trunc(this.length / 2);
        return this.v[index];
    }
    get mode() {
        throw_not_implemented();        
    }
    get min() {
        return Math.min(... this.v);
    }
    get max() {
        return Math.max(... this.v);        
    }
    get best_fit_linear() {
        throw_not_implemented();
    }
    get mse() {
        throw_not_implemented();
    }
    get stddev() {
        throw_not_implemented();        
    }
    fill(filler=0) {
        this.v.fill(filler);
    }
    sort() {
        const vout = this.clone;
        vout.v.sort();
        return vout;
    }
    reverse() {
        const vout = this.clone;
        vout.v.reverse();
        return vout;
    }
    scale(k) {
        const vout = new Vector();
        vout.v = this.v.map(e=>e*k);
        return vout;
    }
    add(other) {
        const vout = this.clone;
        const count = vout.length;
        for(let i=0; i<count; i++) {
            vout.v[i] += other.v[i];
        }
        return vout;
    }
    sub(other) {
        const vout = this.clone;
        const count = vout.length;
        for(let i=0; i<count; i++) {
            vout.v[i] -= other.v[i];
        }
        return vout;        
    }
    mul(other) {
        const vout = this.clone;
        const count = vout.length;
        for(let i=0; i<count; i++) {
            vout.v[i] *= other.v[i];
        }
        return vout;        
    }
    div(other) {
        const vout = this.clone;
        const count = vout.length;
        for(let i=0; i<count; i++) {
            vout.v[i] /= other.v[i];
        }
        return vout;        
    }
    map(callback) {
        const vout = new Vector();
        vout.v = this.v.map(callback);
        return vout;
    }
    reduce(callback) {
        return this.v.reduce(callback);
    }
    iterate(callback) {
        return this.v.forEach(callback);        
    }
    acquire(v) {
        if (v instanceof Vector) v = v.v;
        else if (is_type(v, "string")) v = parseFloat(v);
        else if (is_type(v, "function")) v = v();
        if (!is_iterable(v)) v = [ parseFloat(v) ];
        this.v = [ ... v ];
        return this;
    }
}

// https://en.wikipedia.org/wiki/Complex_number
exports.Complex = class Complex {
    constructor(re, im) {
        this.re = re;
        this.im = im;
    }
    get vec2() {
        return new Vec2(this.re, this.im);
    }
    get polar() {
        return new VecPolar(this.mod, this.arg);
    }
    get mod() {
        return Math.hypot(this.re, this.im);
    }
    get arg() {
        return Math.atan2(this.im, this.re);
    }
    acquire(v) {
        v = new Vec2().acquire(v);
        this.re = v.x;
        this.im = v.im;
        return this;
    }
}

exports.Mtx2x2 = class Mtx2x2 {
    constructor(other=null) {
        this.reset();
        if (other) {
            this.acquire(other);
        }
        this.order = 2;
    }
    get zero() {
        return new Mtx2x2();
    }
    get identity() {
        const m = new Mtx2x2();
        m.m = Mtx2x2.identity();
        return m;
    }
    get clone() {
        const m = new Mtx2x2();
        Matrix.copy(this.m, m.m);
        return m;
    }
    get matrix() {
        const m = new Matrix(2, 2);
        Matrix.copy(this.m, m.m);
        return m;
    }
    reset() {
        this.m = Mtx2x2.zero();
        return this;
    }
    acquire(m) {
        Matrix.copy(m.m, this.m);
        return this;
    }
};

exports.Mtx2x2.zero     = () => [[0,0], [0,0]];
exports.Mtx2x2.identity = () => [[1,0], [0,1]];

exports.Mtx2x3 = class Mtx2x3 {
    constructor(other=null) {
        this.reset();
        if (other) {
            this.acquire(other);
        }
        this.order = 2;
    }
    get zero() {
        return new Mtx2x3();
    }
    get identity() {
        const m = new Mtx2x3();
        m.m = Mtx2x2.identity();
        m.t = Vec2.zero();
        return m;
    }
    get clone() {
        const m = new Mtx2x3();
        Matrix.copy(this.m, m.m);
        Vector.copy(this.t, m.t);
        return m;
    }
    get matrix() {
        const m = new Matrix(2, 3);
        Matrix.copy(this.m, m.m);
        Vector.copy(this.t, m.m[2])
        return m;
    }
    reset() {
        this.m = Mtx2x2.zero();
        this.t = Vec2.zero();
        return this;
    }
    acquire(m) {
        this.reset();
        Matrix.copy(m.m, this.m);
        Vector.copy(m.m[2], this.t);
        return this;
    }
};

exports.Mtx3x3 = class Mtx3x3 {
    constructor(other=null) {
        this.reset();
        if (other) {
            this.acquire(other);
        }
        this.order = 3;
    }
    get zero() {
        return new Mtx3x3();
    }
    get identity() {
        const m = new Mtx3x3();
        m.m = Mtx3x3.identity();
        return m;
    }
    get clone() {
        const m = new Mtx3x3();
        Matrix.copy(this.m, m.m);
        return m;
    }
    get matrix() {
        const m = new Matrix(3, 3);
        Matrix.copy(this.m, m.m);
        return m;
    }
    reset() {
        this.m = Mtx3x3.zero();
        return this;
    }
    acquire(m) {
        this.reset();
        Matrix.copy(m.m, this.m);
        return this;
    }
};

exports.Mtx3x3.zero     = () => [[0,0,0], [0,0,0], [0,0,0]];
exports.Mtx3x3.identity = () => [[1,0,0], [0,1,0], [0,0,1]];

exports.Mtx3x4 = class Mtx3x4 {
    constructor(other=null) {
        this.reset();
        if (other) {
            this.acquire(other);
        }
        this.order = 3;
    }
    get zero() {
        return new Mtx3x4();
    }
    get identity() {
        const m = new Mtx3x4();
        m.m = Mtx3x3.identity();
        m.t = Vec3.zero();
        return m;
    }
    get clone() {
        const m = new Mtx3x4();
        Matrix.copy(this.m, m.m);
        Vector.copy(this.t, m.t);
        return m;
    }
    get matrix() {
        const m = new Matrix(3, 4);
        Matrix.copy(this.m, m.m);
        Vector.copy(this.t, m.m[3])
        return m;
    }
    reset() {
        this.m = Mtx3x3.zero();
        this.t = Vec3.zero();
        return this;
    }
    acquire(m) {
        Matrix.copy(m.m, this.m);
        Vector.copy(m.m[3], this.t);
        return this;
    }
};

exports.Mtx4x4 = class Mtx4x4 {
    constructor(other=null) {
        this.reset();
        if (other) {
            this.acquire(other);
        }
        this.order = 4;
    }
    get zero() {
        return new Mtx4x4();
    }
    get identity() {
        const m = new Mtx4x4();
        m.m = Mtx4x4.identity();
        return m;
    }
    get clone() {
        const m = new Mtx4x4();
        Matrix.copy(this.m, m.m);
        return m;
    }
    get matrix() {
        const m = new Matrix(4, 4);
        Matrix.copy(this.m, m.m);
        return m;
    }
    reset() {
        this.m = Mtx4x4.zero();
        return this;
    }
    acquire(m) {
        Matrix.copy(m.m, this.m);
        return this;
    }
};

exports.Mtx4x4.zero     = () => [[0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0]];
exports.Mtx4x4.identity = () => [[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]];

exports.Matrix = class Matrix {
    constructor(cols=1, rows=1, filler=0) {
        this.m = Matrix.buffer(cols, rows, filler);
        this.order = Math.min(this.rows, this.cols);
    }
    get size() {
        return Matrix.get_size(this.m);
    }
    get rows() {
        return this.m.length || 0;
    }
    get cols() {
        return this.m[0].length || 0;
    }
    get clone() {
        return new Matrix().acquire(this);
    }
    get transpose() {
        const m = new Matrix();
        m.m = Matrix.transpose(this.m);
        return m;
    }
    get inverse() {
        throw_not_implemented();
    }
    get is_square() {
        return this.rows == this.cols;
    }
    get is_row_major() {
        return this.rows > this.cols;
    }
    get is_column_major() {
        return this.rows < this.cols;
    }
    fill(filler=0) {
        Matrix.fill(this.m, 0);
        return this;
    }
    acquire(m) {
        if (Array.isArray(m)) {
            Matrix.copy(m, this.m);
        } else if (Array.isArray(m.m)) {
            Matrix.copy(m.m, this.m);
        } else {
            throw_matrix_size_mismatch();
        }
        return this;
    }
    reset() {
        this.m = Matrix.buffer(this.cols, this.rows);
        return this;
    }
    make_identity() {
        if (!this.is_square) {
            throw_matrix_size_mismatch();
        }
        this.m = Matrix.identity(this.order);
        return this;
    }
    apply(v) {
        const vout = new Vector();
        vout.length = v.length;
        Matrix.apply(this.m, v.v, vout.v);
        return vout;
    }
    concat(m) {
        const mout = new Matrix(m.cols, m.rows);
        Matrix.concat(this.m, m.m, mout.m);
        return mout;
    }
};

exports.Matrix.buffer = function(cols=0, rows=0, filler=0) {
    const m = [];
    while (rows > 0) {
        this.m.push(new Array(cols).fill(filler));
        --rows;
    }
    return m;
}

exports.Matrix.get_size = function(m) {
    const rows = m.length || 0;
    const cols = m[0].length || 0;
    return {rows, cols};
}

exports.Matrix.is_compatible = function(smaller, larger) {
    smaller = Matrix.get_size(smaller);
    larger = Matrix.get_size(larger);
    if (smaller.rows > larger.rows) return false;
    if (smaller.cols > larger.cols) return false;
    return true;
}

exports.Matrix.copy = function(min, mout) {
    if (!Matrix.is_compatible(min, mout)) {
        throw_matrix_size_mismatch();
    }
    const size = Matrix.get_size(min);
    const rows = size.rows;
    const cols = size.cols;
    for (let i=0; i<rows; i++) {
        for (let j=0; j<cols; j++) {
            mout[i][j] = min[i][j];
        }
    }
    return mout;
}

exports.Matrix.fill = function(m, filler) {
    const size = Matrix.get_size(min);
    const rows = size.rows;
    const cols = size.cols;
    for (let i=0; i<rows; i++) {
        for (let j=0; j<cols; j++) {
            m[i][j] = filler;
        }
    }
    return m;
}

exports.Matrix.transpose = function(m) {
    const size = Matrix.get_size(m);
    const rows = size.rows;
    const cols = size.cols;
    const mout = [];
    for (let i=0; i<cols; i++) {
        mout.push(new Array(rows));
    }
    for (let i=0; i<rows; i++) {
        for (let j=0; j<cols; j++) {
            mout[j][i] = m[i][j];
        }
    }
    return mout;
}

exports.Matrix.identity = function(order=1) {
    const m = Matrix.buffer(order, order);
    for (let i=0; i<order; i++) {
        m[i][i] = 1;
    }
    return m;
}

exports.Matrix.apply = function(m, vin, vout) {
    const vinlen = vin.length || 0;
    const voutlen = vout.length || 0;
    const mtxsize = Matrix.get_size(m);
    if (vinlen < mtxsize.cols) throw_vector_size_mismatch();
    if (voutlen < mtxsize.rows) throw_vector_size_mismatch();
    function dot(u, v) {
        let sum = 0;
        for (let i=0; i<u.length; i++) {
            sum += u[i] * v[i];
        }
        return sum;
    }
    for (let i=0; i<mtxsize.rows; i++) {
        vout[i] = dot(vin, m[i]);
    }
    return vout;
}

// Uses premultiply order
exports.Matrix.concat = function(mb, ma, mout) {
    const masize = Matrix.get_size(ma);
    const mbsize = Matrix.get_size(mb);
    const moutsize = Matrix.get_size(mout);
    if (
        (masize.cols < mbsize.rows)   ||
        (masize.rows < moutsize.rows) ||
        (mbsize.cols < moutsize.cols)
    ) {
        throw_matrix_size_mismatch();
    }
    for (let i=0; i<mb.rows; i++) {
        for (let j=0; j<mb.cols; j++) {
            mout[i][j] = 0;
            for (let k=0; k<ma.rows; k++) {
                mout[i][j] += mb[i][k] * ma[k][j];
            }
        }    
    }
    return mout;
}

// Transform for 2-space
exports.Transform2 = class Transform2 {
    // Rotate is scalar
    // Scale and translate are Vec2
    constructor(rotate, scale, translate) {
        this.rotate    = rotate;
        this.scale     = new Vec2(scale.x,     scale.y);
        this.translate = new Vec2(translate.x, translate.y);
    }
    apply(v) {
        return this.mtx3x2.apply(v);
    }
    rotate(v) {
        // TODO...
    }
    scale(v) {
        // TODO...        
    }
    translate(v) {
        // TODO...        
    }
    get mtx2x3() {
        throw_not_implemented();
        /*
        const m = new Mtx2x3();
        return m;
        */
    }
    get mtx3x3() {
        throw_not_implemented();
        /*
        const m = new Mtx3x3();
        return m;
        */
    }
};

// Transform for 3-space
exports.Transform3 = class Transform3 {
    // All three args are Vec3
    constructor(rotate, scale, translate) {
        this.rotate    = new Vec3(rotate.x,    rotate.y,    rotate.z);
        this.scale     = new Vec3(scale.x,     scale.y,     scale.z);
        this.translate = new Vec3(translate.x, translate.y, translate.z);
    }
    apply(v) {
        return this.mtx3x2.apply(v);
    }
    rotate(v) {
        // TODO...
    }
    scale(v) {
        // TODO...        
    }
    translate(v) {
        // TODO...        
    }
    get mtx3x4() {
        throw_not_implemented();
        /*
        const m = new Mtx3x4();
        return m;
        */
    }
    get mtx4x4() {
        throw_not_implemented();
        /*
        const m = new Mtx4x4();
        return m;
        */
    }
};

// This is an generic planar (2D) shape
// For 3D shapes, use the Solid class
exports.Shape = class Shape {
    constructor() {
        this.points = [];
        this.sides = [];
    }
};

// This is like the RECT class from the WinAPI
// https://learn.microsoft.com/en-us/windows/win32/api/windef/ns-windef-rect
exports.Window = class Window {
    constructor(left, top, width, height) {
        this.left = left;
        this.top = top;
        this.width = width;
        this.height = height;
    }
}

// 2D line in slope-intercept form
exports.Line2 = class Line2 {
    constructor(slope=0, intercept=0) {
        this.slope = slope;
        this.intercept = intercept;
    }
    get perpendicular() {
        throw_not_implemented();
    }
    intersect(line) {
        throw_not_implemented();
    }
    // ax+by+c=0
    init_standard(a=0, b=0, c=0) {
        // TODO...
        throw_not_implemented();
    }
    square_distance(pt) {
        throw_not_implemented();
    }
}

// 3D line in standard form
exports.Line3 = class Line3 {
    constructor() {
    }
    get perpendicular() {
        throw_not_implemented();
    }
    intersect(line) {
        throw_not_implemented();
    }
    // ax+by+cz+d=0
    init_standard(a=0, b=0, c=0, d=0) {
        // TODO...
        throw_not_implemented();
    }
    square_distance(pt) {
        throw_not_implemented();
    }
}

// 2D line segment (two points)
exports.LineSeg2 = class LineSeg2 {
    constructor(origin, terminus) {
        this.origin = new Vec2(origin.x, origin.y);
        this.terminus = new Vec2(terminus.x, terminus.y);
    }
    get length() {
        return this.terminus.sub(this.origin).mod;
    }
    get ray() {
        throw_not_implemented();
    }
    get midpoint() {
        throw throw_not_implemented();
    }
    get line() {
        throw_not_implemented();
    }
    // origin=(a, b)
    // terminus=(x, y)
    init_standard(a=0, b=0, x=0, y=0) {
        this.origin = new Vec2(a, b);
        this.terminus = new Vec2(x, y);
    }
}

// 3D line segment (two points)
exports.LineSeg3 = class LineSeg3 {
    constructor(origin, terminus) {
        this.origin = new Vec2(origin.x, origin.y);
        this.terminus = new Vec2(terminus.x, terminus.y);
    }
    get length() {
        return this.terminus.sub(this.origin).mod;
    }
    get ray() {
        throw_not_implemented();
    }
    get midpoint() {
        throw throw_not_implemented();
    }
    get line() {
        throw_not_implemented();
    }
    // origin=(a, b, c)
    // terminus=(x, y, z)
    init_standard(a=0, b=0, c=0, x=0, y=0, z=0) {
        this.origin = new Vec3(a, b, c);
        this.terminus = new Vec3(x, y, z);
    }
}

// https://en.wikipedia.org/wiki/Triangle
exports.Triangle = class Triangle {
    constructor(base, height) {
        this.base = base;
        this.height = height;
    }
    get area() {
        return this.base * this.height * 0.5;
    }    
    get perimeter() {
        throw_not_implemented();
    }
    get midpoint() {
        throw_not_implemented();
    }
}

// https://en.wikipedia.org/wiki/Square
exports.Square = class Square {
    constructor(size) {
        this.size = size;
    }
    get area() {
        return this.size * this.size;
    }    
    get perimeter() {
        return 4 * this.size;
    }
    get midpoint() {
        throw_not_implemented();
    }
}

// https://en.wikipedia.org/wiki/Rectangle
exports.Rectangle = class Rectangle {
    constructor(width, height) {
        this.width = width;
        this.height = height;
    }
    get area() {
        return this.width * this.height;
    }    
    get perimeter() {
        return 2 * (this.width + this.height);
    }
    get midpoint() {
        throw_not_implemented();
    }
}

// Angle is in the lower left corner or upper right corner
// https://en.wikipedia.org/wiki/Parallelogram
exports.Parallelogram = class Parallelogram {
    constructor(width, height, angle) {
        this.width = width;
        this.height = height;
        this.angle = angle;
    }        
    get area() {
        throw_not_implemented()
    }    
    get perimeter() {
        return 2 * (this.width + this.height);
    }
    get midpoint() {
        throw_not_implemented();
    }
}

// https://en.wikipedia.org/wiki/Rhombus
exports.Rhombus = class Rhombus {
    // Args are points
    constructor(a, b, c, d) {
        this.points = this.acquire([
            a, b, c, d
        ]);
    }
    get a() {
        return this.points[0];
    }
    get b() {
        return this.points[1];        
    }
    get c() {
        return this.points[2];        
    }
    get d() {
        return this.points[3];        
    }
    get area() {
        throw_not_implemented()
    }    
    get perimeter() {
        return this.top + this.bottom + this.left + this.right;
    }
    get midpoint() {
        throw_not_implemented();
    }
    sort() {
        // TODO: this needs to pick a point to serve as
        // the lower left corner, then walk CCW around
        // the perimeter assigning points in winding order
        const mod2pi = radians => {
            return wrap(radians, Math.PI*2);
        }
        const cmp = (a, b) => {
            const aa = mod2pi(a.arg);
            const ba = mod2pi(b.arg);
            return a - b;
        }
        this.points = this.points.sort(cmp);
        return this;
    }
    acquire(points) {
        while (points.length < 4) points.push(new Vec2());
        this.points = points.slice(0, 4).map(pt=>pt.clone);
        return this.sort();
    }
}

// https://en.wikipedia.org/wiki/Conic_section
exports.Circle = class Circle {
    constructor(radius) {
        this.radius = radius;
    }
    get area() {
        return this.radius*this.radius*Math.PI;
    }
    get diameter() {
        return this.radius*2;
    }
    get circumference() {
        return this.diameter*Math.PI;
    }
    get midpoint() {
        throw_not_implemented();
    }
}

// https://en.wikipedia.org/wiki/Conic_section
exports.Ellipse = class Ellipse {
    constructor(radius, aspect) {
        this.radius = radius;
        this.aspect = aspect;
    }                
    get area() {
        throw_not_implemented();
    }
    get diameter() {
        throw_not_implemented();
    }
    get circumference() {
        throw_not_implemented();
    }
    get midpoint() {
        throw_not_implemented();
    }
}

// https://en.wikipedia.org/wiki/Parabola
// https://en.wikipedia.org/wiki/Conic_section
// https://byjus.com/maths/standard-equations-of-parabola/
exports.Parabola = class Parabola {
    constructor(vertex, focus) {
        this.vertex = new Vec2(vertex.x, vertex.y);
        this.focus = new Vec2(focus.x, focus.y);
    }
}

// https://en.wikipedia.org/wiki/Conic_section
// https://www.cuemath.com/geometry/hyperbola/
exports.Hyperbola = class Hyperbola {
    constructor(vertex, focus) {
        this.vertex = new Vec2(vertex.x, vertex.y);
        this.focus = new Vec2(focus.x, focus.y);
    }
}

// This is an generic polyhedral (3D) shape
// For 2D shapes, use the Shape class
// 
exports.Solid = class Solid {
    constructor() {
        this.points = [];
        this.sides = [];
    }
};

// This is a ray in 2D space
// https://en.wikipedia.org/wiki/Line_(geometry)#Ray
exports.Ray2 = class Ray2 {
    constructor(origin, direction) {
        this.origin = new Vec2(origin.x, origin.y, origin.z);
        this.direction = new Vec2(direction.x, direction.y, direction.z);
    }
};

// This is a ray in 3D space
// https://en.wikipedia.org/wiki/Line_(geometry)#Ray
// https://en.wikipedia.org/wiki/Ray_(optics)
exports.Ray3 = class Ray3 {
    constructor(origin, direction) {
        this.origin = new Vec3(origin.x, origin.y, origin.z);
        this.direction = new Vec3(direction.x, direction.y, direction.z);
    }
};

// https://en.wikipedia.org/wiki/Euclidean_plane
exports.Plane = class Plane {
    constructor(normal, distance) {
        this.normal = new Vec3(normal.x, normal.y, normal.z);
        this.distance = distance;
    }
    // ax+by+cz+d=0
    init_standard(a, b, c, d) {
        // TODO...
    }
};

// https://en.wikipedia.org/wiki/Cube
exports.Cube = class Cube {
    constructor(size) {
        this.size = size;
    }
};

// https://en.wikipedia.org/wiki/Cuboid
exports.Cuboid = class Cuboid {
    constructor(width, height, depth) {
        this.width = width;
        this.height = height;
        this.depth = depth;
    }
};

// https://en.wikipedia.org/wiki/Sphere
exports.Sphere = class Sphere {
    constructor(radius) {
        this.radius = radius;
    }
};

// https://en.wikipedia.org/wiki/Spheroid
exports.Spheroid = class Spheriod {
    constructor(radius, aspect) {
        this.radius = radius;
        this.aspect = aspect;
    }
};

// https://en.wikipedia.org/wiki/Ellipsoid
exports.Ellipsoid = class Ellipsoid {
    constructor(xradius, yradius, zradius) {
        this.xradius = xradius;
        this.yradius = yradius;
        this.zradius = zradius;
    }
};

// This is a regular polyhedron
// https://en.wikipedia.org/wiki/Polyhedron
// https://www.cuemath.com/geometry/polyhedron/
exports.Polyhedron = class Polyhedron {
    constructor(radius, nsides) {
        this.radius = radius;
        this.nsides = nsides;
    }
};

// Represents a polygon in 3D
exports.PolygonPatch = class PolygonPatch {
    constructor() {
        throw_not_implemented();
    }
};

// Represents a triangle in 3D
exports.TrianglePatch = class TrianglePatch {
    constructor() {
        throw_not_implemented();
    }
};

// Represents a rectangle in 3D
exports.RectanglePatch = class RectanglePatch {
    constructor() {
        throw_not_implemented();
    }
};

// This is a cylinder-like solid, but faceted
exports.CylinderPost = class CylinderPost {
    constructor() {
        throw_not_implemented();
    }
};

// https://en.wikipedia.org/wiki/Cylinder
exports.Cylinder = class Cylinder {
    constructor() {
        throw_not_implemented();
    }
};

// https://en.wikipedia.org/wiki/Cone
exports.Cone = class Cone {
    constructor() {
        throw_not_implemented();
    }
}

// https://en.wikipedia.org/wiki/Pyramid_(geometry)
exports.Pyramid = class Pyramid {
    constructor() {
        throw_not_implemented();
    }
};

// Planar (flat) disc shape with optional
// hole in center
// https://en.wikipedia.org/wiki/Disk_(mathematics)
// https://en.wikipedia.org/wiki/Hard_disk_drive_platter
exports.Platter = class Platter {
    constructor() {
        throw_not_implemented();
    }
};

// https://en.wikipedia.org/wiki/Toroid
exports.Torus = class Torus {
    constructor() {
        throw_not_implemented();
    }
};

// https://en.wikipedia.org/wiki/Toroid
exports.Toroid = class Toroid {
    constructor() {
        throw_not_implemented();
    }
};

// https://en.wikipedia.org/wiki/Platonic_solid
// https://chat.openai.com/c/53542719-aa4e-4b26-903c-a1b2d5901fcb
exports.PlatonicSolid = class PlatonicSolid {
    constructor(faces=4) {
        this.faces = faces;
    }
    get types() {
        const type = (title, faces, builder) => { title, faces };
        return [
            type("tetrahedron", 4),
            type("cube", 6),
            type("octahedron", 8),
            type("dodecahedron", 12),
            type("icosahedron", 20)
        ];
    }
    build() {
        switch(this.faces) {
        case 4: return this.build_tetra();
        case 6: return this.build_cube();
        case 8: return this.build_octa();
        case 12: return this.build_dodeca();
        case 20: return this.build_icosa();
        default:
            throw new Error("There is no platonic solid with this number of faces.");
        }
    }
    build_tetra() {
        throw_not_implemented();
    }
    build_cube() {
        throw_not_implemented();
    }
    build_octa() {
        throw_not_implemented();
    }
    build_dodeca() {
        throw_not_implemented();
    }
    build_icosa() {
        throw_not_implemented();
    }
};

// https://www.kristakingmath.com/blog/direction-cosines-direction-angles
exports.DirectionAngles = class DirectionAngles {
    constructor(other=null) {
        this.acquire(other);
    }
    get cosines() {
        const n = this.length;        
        const t = new Array(n);
        for (let i=0; i<n; i++) {
            t[i] = Math.cos(this.angles[i]);
        }
        return t;
    }
    set cosines(source) {
        const n = source.length;        
        this.angles = new Array(n);
        for (let i=0; i<n; i++) {
            this.angles[i] = Math.acos(source[i]);
        }
        return this;
    }
    get clone() {
        return new DirectionAngles(this);
    }
    get length() {
        return this.angles.length;
    }
    set length(n) {
        this.angles = new Array(n);
        return this;
    }
    fill(filler=0) {
        this.angles = new Array(this.length).fill(filler);
    }
    project(distance) {
        const t = this.cosines;
        const n = this.length;
        const v = new Vector().length = n;
        for (let i=0; i<n; i++) {
            v.v[i] = t[i] * distance;
        }
        return v;
    }
    acquire(other=null) {
        if (other instanceof DirectionAngles) {
            other = other.angles;
        } else if (is_type(other, 'string')) {
            other = [ parseFloat(other) ];
        } else if (is_type(other, 'function')) {
            other = other();
        } else if (is_iterable(other)) {
            other = [ ... other ];
        } else if (!isNaN(other)) {
            other = [ parseFloat(other) ];
        } else {
            other = [];
        }
        this.angles = [ ... other ];
        return this;
    }
};

