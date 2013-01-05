fn main() {
	let n = int::from_str(os::args()[1]).get();

	let bodies = n_body_system::new();
	io::println(fmt!("%.9f", bodies.energy() as float));
	let mut i = 0;
	while i < n {
		bodies.advance(0.01f64);
		i += 1;
	}
	io::println(fmt!("%.9f", bodies.energy() as float));
}

struct n_body_system {
	priv bodies: ~[@body]
}

impl n_body_system {
	static fn new() -> n_body_system {
		let result = n_body_system { bodies: ~[sun(), jupiter(), saturn(), uranus(), neptune()] };

		let mut px = 0f64;
		let mut py = 0f64;
		let mut pz = 0f64;
		let mut i = 0;
		while i < result.bodies.len() {
			px += result.bodies[i].vx * result.bodies[i].mass;
			py += result.bodies[i].vy * result.bodies[i].mass;
			pz += result.bodies[i].vz * result.bodies[i].mass;
			i += 1;
		}
		result.bodies[0].offsetMomentum(px, py, pz);

		result
	}

	fn advance(&self, dt: f64) {
		let mut i = 0;
		while i < self.bodies.len() {
			let ibody = self.bodies[i];
			let mut j = i + 1;
			while j < self.bodies.len() {
				let jbody = self.bodies[j];
				let dx = ibody.x - jbody.x;
				let dy = ibody.y - jbody.y;
				let dz = ibody.z - jbody.z;

				let dsquared = dx * dx + dy * dy + dz * dz;
				let distance = f64::sqrt(dsquared);
				let mag = dt / (dsquared * distance);

				ibody.vx -= dx * jbody.mass * mag;
				ibody.vy -= dy * jbody.mass * mag;
				ibody.vz -= dz * jbody.mass * mag;

				jbody.vx += dx * ibody.mass * mag;
				jbody.vy += dy * ibody.mass * mag;
				jbody.vz += dz * ibody.mass * mag;

				j += 1;
			}
			i += 1;
		}

		let mut n = 0;
		while n < self.bodies.len() {
			let body = self.bodies[n];
			body.x += dt * body.vx;
			body.y += dt * body.vy;
			body.z += dt * body.vz;

			n += 1;
		}
	}

	fn energy(&self) -> f64 {
		let mut e = 0f64;

		let mut i = 0;
		while i < self.bodies.len() {
			let ibody = self.bodies[i];
			e += 0.5f64 * ibody.mass * (ibody.vx * ibody.vx + ibody.vy * ibody.vy + ibody.vz * ibody.vz);
			let mut j = i + 1;
			while j < self.bodies.len() {
				let jbody = self.bodies[j];
				let dx = ibody.x - jbody.x;
				let dy = ibody.y - jbody.y;
				let dz = ibody.z - jbody.z;

				let distance = f64::sqrt(dx * dx + dy * dy + dz * dz);
				e -= (ibody.mass * jbody.mass) / distance;
				j += 1;
			}
			i += 1;
		}
		e
	}
}

const pi: f64 = 3.141592653589793f64;
const solar_mass: f64 = 4f64 * pi * pi;
const days_per_year: f64 = 365.24f64;

struct body {
	mut x: f64,
	mut y: f64,
	mut z: f64,
	mut vx: f64,
	mut vy: f64,
	mut vz: f64,
	mass: f64
}

impl body {
	static fn new(x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64, mass: f64) -> body {
		body { x: x, y: y, z: z, vx: vx, vy: vy, vz: vz, mass: mass }
	}

	fn offsetMomentum(&self, px: f64, py: f64, pz: f64) {
		self.vx = -px / solar_mass;
		self.vy = -py / solar_mass;
		self.vz = -pz / solar_mass;
	}
}

fn jupiter() -> @body {
	@body::new(
		4.84143144246472090e+00f64,
		-1.16032004402742839e+00f64,
		-1.03622044471123109e-01f64,
		1.66007664274403694e-03f64 * days_per_year,
		7.69901118419740425e-03f64 * days_per_year,
		-6.90460016972063023e-05f64 * days_per_year,
		9.54791938424326609e-04f64 * solar_mass)
}

fn saturn() -> @body {
	@body::new(
		8.34336671824457987e+00f64,
		4.12479856412430479e+00f64,
		-4.03523417114321381e-01f64,
		-2.76742510726862411e-03f64 * days_per_year,
		4.99852801234917238e-03f64 * days_per_year,
		2.30417297573763929e-05f64 * days_per_year,
		2.85885980666130812e-04f64 * solar_mass)
}

fn uranus() -> @body {
	@body::new(
		1.28943695621391310e+01f64,
		-1.51111514016986312e+01f64,
		-2.23307578892655734e-01f64,
		2.96460137564761618e-03f64 * days_per_year,
		2.37847173959480950e-03f64 * days_per_year,
		-2.96589568540237556e-05f64 * days_per_year,
		4.36624404335156298e-05f64 * solar_mass)
}

fn neptune() -> @body {
	@body::new(
		1.53796971148509165e+01f64,
		-2.59193146099879641e+01f64,
		1.79258772950371181e-01f64,
		2.68067772490389322e-03f64 * days_per_year,
		1.62824170038242295e-03f64 * days_per_year,
		-9.51592254519715870e-05f64 * days_per_year,
		5.15138902046611451e-05f64 * solar_mass)
}

fn sun() -> @body {
	@body::new(
		0f64,
		0f64,
		0f64,
		0f64,
		0f64,
		0f64,
		solar_mass)
}
