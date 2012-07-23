
fn main(args: ~[str]) {

	let ALU: str = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGG\
						 GAGGCCGAGGCGGGCGGATCACCTGAGGTCAGGAGTTCGAGA\
						 CCAGCCTGGCCAACATGGTGAAACCCCGTCTCTACTAAAAAT\
						 ACAAAAATTAGCCGGGCGTGGTGGCGCGCGCCTGTAATCCCA\
						 GCTACTCGGGAGGCTGAGGCAGGAGAATCGCTTGAACCCGGG\
						 AGGCGGAGGTTGCAGTGAGCCGAGATCGCGCCACTGCACTCC\
						 AGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAA";

	let IUB: ~float_prob_freq = ~float_prob_freq(
				@['a',  'c',  'g',  't',
	              'B',  'D',  'H',  'K',
	              'M',  'N',  'R',  'S',
	              'V',  'W',  'Y'],
	        	~[0.27f32, 0.12f32, 0.12f32, 0.27f32,
	              0.02f32, 0.02f32, 0.02f32, 0.02f32,
	              0.02f32, 0.02f32, 0.02f32, 0.02f32,
	              0.02f32, 0.02f32, 0.02f32], 42);

	let n = if args.len() > 1 { uint::from_str(args[1]).get() } else { 1000u };

	let out = io::stdout();
	make_repeat_fasta("ONE", "Homo sapiens alu", ALU, n * 2, out);
	make_random_fasta("TWO", "IUB ambiguity codes", IUB, n * 3, out);

	let homo_sapiens: ~float_prob_freq = ~float_prob_freq(
				@['a', 'c', 'g', 't'],
				~[0.3029549426680f32,
	              0.1979883004921f32,
	              0.1975473066391f32,
	              0.3015094502008f32], IUB.last);
	make_random_fasta("THREE", "Homo sapiens frequency", homo_sapiens, n * 5, out);
	out.flush();
}

const IM: int = 139968;
const IA: int = 3877;
const IC: int = 29573;

const line_length: uint = 60;
const buffer_size: uint = (line_length + 1) * 1024;

fn make_random_fasta(id: str, desc: str, fpf: ~float_prob_freq, n_chars: uint, writer: io::writer) {
	let buffer = vec::to_mut(vec::from_elem(buffer_size, 0u8));
	if buffer.len() % (line_length + 1) != 0 {
		log(error, "buffer size must be a multiple of line length (including line break)");
		fail;
	}

	let desc_str = #fmt(">%s %s\n", id, desc);
	writer.write(str::bytes(desc_str));

	let mut buffer_index = 0;
	let mut n = n_chars;
	while n > 0 {
		let chunk_size = if n >= line_length { line_length } else { n };
		if buffer_index == buffer_size {
			writer.write(buffer);
			buffer_index = 0;
		}
		buffer_index = fpf.select_random_into_buffer(buffer, buffer_index, chunk_size);
		buffer[buffer_index] = '\n' as u8;
		buffer_index += 1;
		n -= chunk_size;
	}

	writer.write(vec::view(buffer, 0, buffer_index));
}

fn make_repeat_fasta(id: str, desc: str, alu: str, n_chars: uint, writer: io::writer) {
	let alu_bytes = str::bytes(alu);
	let mut alu_index = 0;

	let buffer = vec::to_mut(vec::from_elem(buffer_size, 0u8));
	if buffer.len() % (line_length + 1) != 0 {
		log(error, "buffer size must be a multiple of line length (including line break)");
		fail;
	}

	let desc_str = #fmt(">%s %s\n", id, desc);
	writer.write(str::bytes(desc_str));

	let mut buffer_index = 0;
	let mut n = n_chars;
	while n > 0 {
		let chunk_size = if n >= line_length { line_length } else { n };
		if buffer_index == buffer_size {
			writer.write(buffer);
			buffer_index = 0;
		}

		let mut i = 0;
		while i < chunk_size {
			if alu_index == alu_bytes.len() {
				alu_index = 0;
			}
			buffer[buffer_index] = alu_bytes[alu_index];
			buffer_index += 1;
			alu_index += 1;
			i += 1;
		}
		buffer[buffer_index] = '\n' as u8;
		buffer_index += 1;

		n -= chunk_size;
	}

	writer.write(vec::view(buffer, 0, buffer_index));
}

class float_prob_freq {
	let mut last: int;
	priv {
		let chars: @[char];
		let probs: ~[mut f32];
	}

	new(chars: @[char], probs: ~[f32], last: int) {
		self.last = last;
		self.chars = chars;
		self.probs = vec::to_mut(copy probs);

		self.make_cumulative();
	}

	fn make_cumulative() {
		let mut cp = 0.0f64;
		let mut i = 0;
		while i < self.probs.len() {
			cp += self.probs[i] as f64;
			self.probs[i] = cp as f32;
			i += 1;
		}
	}

	fn select_random_into_buffer(buffer: ~[mut u8], buffer_index: uint, n_random: uint) -> uint {
		let len = self.probs.len();

		let mut bi = buffer_index;
		let mut r_index = 0;
		while r_index < n_random {
			let r = self.random(1f32);
			let mut i = 0;
			let mut b = false;
			while i < len {
				if r < self.probs[i] {
					buffer[bi] = self.chars[i] as u8;
					bi += 1;
					b = true;
					break;
				}
				i += 1;
			}
			if !b {
				buffer[bi] = self.chars[len - 1] as u8;
				bi += 1;
			}
			r_index += 1;
		}

		bi
	}

	fn random(max: f32) -> f32 {
		let one_over_IM = 1f32 / (IM as f32);
		self.last = (self.last * IA + IC) % IM;
		max * (self.last as f32) * one_over_IM
	}
}
