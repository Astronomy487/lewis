let message_box = document.querySelector("#message_box");
let title_box = document.querySelector("#title_box");
function parse_and_run(formula, charge) {
	while (document.querySelector("main")) document.querySelector("main").remove();
	title_box.innerHTML = "";
	message_box.innerText = "computing lewis structure...";
	setTimeout(async function(){
		if (formula == "") return;
		//umm maybe sanitize input even more. like what if we are given something that isnt even a formula
		let results = solve(formula, charge);
		for (let i = 0; i < 10000; i++) {
			let new_results = solve(formula, charge);
			if (new_results.problematic_score < results.problematic_score) results = new_results;
		}
		if (results.problematic_score >= 1000000-1) {
			message_box.innerText = "couldn't compute a good lewis structure\nproblems:\n"+results.significant_problems.join("\n");
		} else {
			message_box.innerText = "drawing lewis structure...";
			//do the draw thing
			setTimeout(async function(){
				let attempts = 0;
				let best = {overlap: Number.MAX_SAFE_INTEGER};
				//for (let attempts = 0; attempts < 100 || best.overlap > 100000; attempts++) {
				while (best.overlap > 100000 && attempts < 10000) {
					let new_attempt = try_draw(results.molecule);
					if (new_attempt.overlap < best.overlap) best = new_attempt;
					attempts++;
				}
				if (best.overlap > 100000) {
					message_box.innerText = "couldn't draw a good lewis structure.";
				} else {
					message_box.innerText = "";
					document.body.appendChild(best.main);
					title_box.innerHTML = make_formula_pretty(formula) + (charge==0 ? "" : " <sup>"+make_charge_pretty(charge)+"</sup>") + "<br>" + symbol_to_molar_mass(formula).toFixed(1) + " g/mol";
				}
			});
		}
	});
}

function solve(formula, charge) {
	let atom_counts = symbol_to_atom_counts(formula);
	let molecule = {
		atoms: [],
		net_formal_charge: charge
	};
	for (let element of Object.keys(atom_counts)) for (let i = 0; i < atom_counts[element]; i++) {
		molecule.atoms.push({symbol: element, bonds: []});
	}
	
	//find how many electrons each atom brings
	for (let atom of molecule.atoms) {
		atom.atom_object = symbol_to_element_object(atom.symbol);
		let group = atom.atom_object.group;
		let brings = 2; //for the d and f blocks. very awkward
		if (group <= 2) brings = group;
		if (group >= 13) brings = group - 10;
		atom.brings = brings;
	}

	for (let atom of molecule.atoms) atom.lone_electrons = 0;
	
	//sort the molecules so that carbon-core things are in the middle
	molecule.atoms.sort((a, b) => Math.random()-0.5);
	//molecule.atoms.sort((a, b) => heuristic_centrality(b) - heuristic_centrality(a));
	for (let i = 1; i < molecule.atoms.length; i++) molecule.atoms[i].bonds.push(Math.floor(Math.random()*i/1)); //THIS IS BAD!!! the way we connect this molecule shouldnt just be random you idiot
	
	/*
	idea: instead of flattening a formula as we recieve it,
	make a LIST of atoms that are present
	the order roughly comes from the formula
	you can create a queue of sorts of like The Last Atoms That Were Just Mentioned
	and in general one of the last ones should be the one this new atom connects to
	The Order of Formulae tend to mean things!!!
	maybe this will enable large hydrocarbon things to Work
	*/

	//turn bonds from index arrays into actual {atoms: [a, b], order: 0}
	//the two sides of the atom most hold the same bond object!!!
	//atoms list within bond object is unordered. theres a helper function to find the Other atom if u give it the atom you are currently on
	for (let i = 0; i < molecule.atoms.length; i++) {
		for (let j = i+1; j < molecule.atoms.length; j++) {
			if (molecule.atoms[i].bonds.includes(j) || molecule.atoms[j].bonds.includes(i)) {
				//i and j ARE connected. sever their fake number link
				molecule.atoms[i].bonds = remove_value(molecule.atoms[i].bonds, j);
				molecule.atoms[j].bonds = remove_value(molecule.atoms[j].bonds, i);
				//make the real links
				let the_bond = {
					atoms: [molecule.atoms[i], molecule.atoms[j]],
					order: 0
				};
				molecule.atoms[i].bonds.push(the_bond);
				molecule.atoms[j].bonds.push(the_bond);
			}
		}
	}

	//sort the real actual atom list for number of bonds (roughly corresponds to centrality. goes outer to inner)
	molecule.atoms.sort((a, b) => a.bonds.length - b.bonds.length);

	//find number of electrons that need to be assigned at some point
	let electrons_in_molecule = 0 - molecule.net_formal_charge;
	for (let atom of molecule.atoms) electrons_in_molecule += atom.brings;
	let electrons_to_assign = electrons_in_molecule;

	//put electrons into bonds to get them to order one
	for (let atom of molecule.atoms) for (let bond of atom.bonds) {
		if (bond.order == 0) {
			bond.order++;
			electrons_to_assign -= 2;
		}
	}

	//going from outside to inside, just assign electrons to satisfy octets
	for (let atom of molecule.atoms) {
		while (!satisfied_octet(atom) && electrons_to_assign) {
			atom.lone_electrons++;
			electrons_to_assign--;
		}
	}

	//if we have leftover electrons, try really hard to find expandable octets. especially among positive formal charge atoms
	while (electrons_to_assign >= 2) {
		let atoms_from_positive_to_negative_charge = molecule.atoms.toSorted((b, a) => atom_formal_charge(a) - atom_formal_charge(b));
		//we have this list. we should add 2 to the most positive atom, if it can take it.
		//if it can't, we keep going through list until that can happen.
		//if we ever go through whole list without making any improvements, break out of this while
		let any_improvements = false;
		for (let atom of atoms_from_positive_to_negative_charge) {
			if (how_much_can_octet_expand(atom) >= 2 && electrons_to_assign >= 2) {
				atom.lone_electrons += 2;
				electrons_to_assign -= 2;
				any_improvements = true;
				break;
			}
		}
		if (!any_improvements) break;
	}
		
	//sweep around atoms. if i have negative charge, >=2 lone electrons, and know someone with positive charge, let's turn that into a stronger bond. (unless that other atom cannot have more)
	{
		let any_effect = true;
		while (any_effect) {
			any_effect = false;
			for (let atom of molecule.atoms) {
				for (let bond of atom.bonds) {
					let partner = bond_other_atom(atom, bond);
					if (atom_formal_charge(partner) > 0 && atom_formal_charge(atom) < 0 && atom.lone_electrons >= 2 && can_take_more_electrons(partner)>=2) {
						atom.lone_electrons -= 2;
						bond.order++;
						any_effect = true;
					}
				}
			}
		}
	}
	
	//OK FINAL CHECK: if anything sucks super bad. um sorry, try again i guess
	let problematic_score = 0;
	let significant_problems = [];
	//unassigned electron
	if (electrons_to_assign != 0) {
		problematic_score += 1000000;
		significant_problems.push("we literally just can't assign the number of electrons");
	}
	for (let atom of molecule.atoms) {
		//literally exceeding octet
		if (too_many_electrons(atom)) {
			problematic_score += 1000000;
			significant_problems.push(atom.symbol + " exceeds octet");
		}
		//undersatisfied octet
		if (!satisfied_octet(atom)) {
			problematic_score += 1000000;
			significant_problems.push(atom.symbol + " has unsatisfied octet");
		}
		//formal charges (tiny problem)
		problematic_score += Math.abs(atom_formal_charge(atom)) * 10;
		//using expanded octets (even tinier problem)
		problematic_score += Math.max(atom_electrons_felt(atom)-8, 0);
	}
	//TODO else check that negative charges are on the electronegative atoms
	
	return {molecule: molecule, problematic_score: problematic_score, significant_problems: significant_problems};
}

function make_charge_pretty(formal_charge) {
	if (formal_charge == 0) return "";
	return (Math.abs(formal_charge)==1 ? "" : Math.abs(formal_charge)) + (formal_charge>0 ? "+" : "–")
}

function try_draw(molecule) {
	let main = document.createElement("main");
	let first_choice_atom = molecule.atoms[0]; //this makes the first choice an extraneous atom (like a hydrogen or halogen) I did this because it makes the animation look cool
	first_choice_atom.x = 0;
	first_choice_atom.y = 0;
	let draw_queue = [{
		atom: first_choice_atom,
		//coming_from_angle: Math.random() * 2 * Math.PI,
		coming_from_angle: Math.random()<0.5 ? 0 : Math.PI,
		previous_atom: null
	}];
	//function draw_atom(atom, coming_from_angle, x, y, previous_atom) {
	let object_number = 0; //used to make ripply animation :)
	while (draw_queue.length) {
		object_number += 2 / molecule.atoms.length;
		//STUPID QUEUE SYSTEM BECAUSE MY RECURSIVE FUNCTIONS ARENT WORKING
		let queue_item = draw_queue.shift();
		let [atom, coming_from_angle, x, y, previous_atom] = [queue_item.atom, queue_item.coming_from_angle, queue_item.atom.x, queue_item.atom.y, queue_item.previous_atom];
		//make an element
		atom.element = main.appendChild(document.createElement("div"));
		atom.element.setAttribute("class", "atom");
		atom.element.innerText = atom.symbol;
		//if (atom.atom_object.color) atom.element.style.outline = "solid 0.25rem " + atom.atom_object.color;
		if (atom.atom_object.color) atom.element.style.backgroundColor = atom.atom_object.color;
		let formal_charge = atom_formal_charge(atom);
		if (formal_charge != 0) {
			let charge_symbol = main.appendChild(document.createElement("div"));
			charge_symbol.setAttribute("class", "charge");
			charge_symbol.style.left = x+1+"rem";
			charge_symbol.style.top = y-0.5+"rem";
			charge_symbol.innerText = make_charge_pretty(formal_charge);
		}
		atom.element.style.left = x + "rem";
		atom.element.style.top = y + "rem";
		atom.element.style.animationDuration = object_number+"s";
		//atom.element.style.backgroundColor = atom.atom_object.color;
		let steric = steric_number(atom);
		//compile list of other things to draw
		things_to_draw = atom.bonds.filter((bond) => !bond.atoms.includes(previous_atom));
		shuffle(things_to_draw);
		for (let i = 0; i < atom.lone_electrons - 1; i += 2) things_to_draw.push(2);
		if (atom.lone_electrons % 2) things_to_draw.push(1);
		//create the set of new angles where are things to draw should be put
		let angles_to_do = [];
		for (let i = 0; i < things_to_draw.length; i++) angles_to_do.push((coming_from_angle + Math.PI + Math.PI*2*(i+1)/steric) % (Math.PI * 2));
		shuffle(angles_to_do);
		for (let thing_to_draw of things_to_draw) {
			let distance_amplitude = typeof(thing_to_draw) == "number" ? 1.5 : 5;
			//find the best angle to draw this thing at
			let angle = 0;
			let largest_sum_of_distances_to_others = 0;
			for (let possible_angle of angles_to_do) {
				let [possible_x, possible_y] = [x + distance_amplitude * Math.cos(possible_angle), y + distance_amplitude * Math.sin(possible_angle)]; //position of this link
				let sum_of_distances_to_others = 0;
				for (let other_atom of molecule.atoms) if (other_atom.x != undefined) sum_of_distances_to_others += Math.pow(Math.hypot(other_atom.x - possible_x, other_atom.y - possible_y), 2); //square distances because it feels right
				if (sum_of_distances_to_others > largest_sum_of_distances_to_others) [angle, largest_sum_of_distances_to_others] = [possible_angle, sum_of_distances_to_others];
			}
			angles_to_do = remove_value(angles_to_do, angle);
			if (typeof(thing_to_draw) == "number") {
				//draw a lone pair/electron here
				let lone = main.appendChild(document.createElement("div"));
				lone.setAttribute("class", "lone");
				lone.innerText = {2:":",1:"·"}[thing_to_draw];
				lone.style.left = x + distance_amplitude * Math.cos(angle) + "rem";
				lone.style.top = y + distance_amplitude * Math.sin(angle) + "rem";
				lone.style.rotate = angle + "rad";
				lone.style.animationDuration = object_number+"s";
			} else {
				//draw line and the next atom
				let new_x = x + distance_amplitude * Math.cos(angle);
				let new_y = y + distance_amplitude * Math.sin(angle);
				let order = thing_to_draw.order;
				for (let o = 0; o < order; o++) {
					let altitude = (o - (order-1)*0.5) * 0.5;
					let perp_angle = angle + Math.PI / 2;
					//all the coords will be moved by alt * cis(perp_angle)
					make_hr(x + altitude * Math.cos(perp_angle), y + altitude * Math.sin(perp_angle), new_x + altitude * Math.cos(perp_angle), new_y + altitude * Math.sin(perp_angle));
				}
				let partner = bond_other_atom(atom, thing_to_draw);
				//WHEN WE PUT THINGS IN DRAW QUEUE!!! make their atom objects have their x position so we can start using that info already!
				partner.x = new_x;
				partner.y = new_y;
				draw_queue.push({
					atom: partner,
					coming_from_angle: angle,
					previous_atom: atom
				});
			}
		}
	}
	function make_hr(x1, y1, x2, y2, ratio) {
		let centerx = (x1+x2) / 2;
		let centery = (y1+y2) / 2;
		let hr = main.appendChild(document.createElement("hr"));
		hr.style.top = centery+0.45 + "rem";
		hr.style.left = centerx-1.55 + "rem";
		hr.style.width = Math.hypot(x2-x1, y2-y1) + "rem";
		//this altitude is messed up because of transformation orders. i hate that for me
		//hr.style.transform = "translateY("+altitude+"rem) rotate("+(Math.PI/2-Math.atan2(x2-x1, y2-y1))+"rad)";
		hr.style.rotate = Math.PI/2 - Math.atan2(x2-x1, y2-y1) + "rad";
		//hr.style.borderWidth = order*order * 0.0625 + "rem";
		hr.style.animationDuration = object_number+"s";
	}
	let overlap = 0;
	for (let i = 0; i < molecule.atoms.length; i++) for (let j = i+1; j < molecule.atoms.length; j++) {
		let [atom_i, atom_j] = [molecule.atoms[i], molecule.atoms[j]];
		let distance = Math.hypot(atom_i.x - atom_j.x, atom_i.y - atom_j.y);
		if (distance < 3) overlap += 100000;
		if (distance < 5) overlap += 5;
		if (distance < 7) overlap++;
		if (distance < 9) overlap++;
	}
	//make height-width of main to be 2rem + the difference between most and least coordinate
	let [least_x, most_x, least_y, most_y] = [0, 0, 0, 0];
	for (let atom of molecule.atoms) {
		least_x = Math.min(least_x, atom.x);
		most_x = Math.max(most_x, atom.x);
		least_y = Math.min(least_y, atom.y);
		most_y = Math.max(most_y, atom.y);
	}
	main.style.transform = "translate(-50%, -50%) translate("+(-1 - 0.5*most_x - 0.5*least_x)+"rem, "+(-1 - 0.5*most_y - 0.5*least_y)+"rem)";
	return {main: main, overlap: overlap};
}

function shuffle(r){for(var a=r.length-1;0<a;a--){var f=Math.floor(Math.random()*(a+1)),o=r[a];r[a]=r[f],r[f]=o}}

function satisfied_octet(atom) {
	let electrons_felt = atom.lone_electrons;
	for (let bond of atom.bonds) electrons_felt += bond.order * 2;
	if (atom.atom_object.period == 1) return electrons_felt == 2;
	if (atom.atom_object.group == 13 && electrons_felt == 6) return true;
	if (atom.atom_object.group == 12 && electrons_felt == 4) return true;
	if (atom.atom_object.z <= 14) return electrons_felt == 8;
	return electrons_felt >= 8; //is there any upper limit on the expanded octet rule???
}

function how_much_can_octet_expand(atom) {
	if (atom.atom_object.z <= 14) return 0;
	let electrons_felt = atom.lone_electrons;
	for (let bond of atom.bonds) electrons_felt += bond.order * 2;
	return 18 - electrons_felt;
}

function heuristic_centrality(atom) { //based off of BRINGS. ranks fluorine(1) = hydrogen(1) < oxygen(2) < nitrogen < carbon. doesn't consider expanded octets
	let bigness_component = atom.atom_object.z * atom.atom_object.z * 0.00005;
	if (atom.atom_object.period == 1) return 2 - atom.brings + bigness_component;
	return 8 - atom.brings + bigness_component;
	//maybe for super big atoms, up the number?? we want Xe to seem more central in XeI2
}

function can_take_more_electrons(atom) { //returns number of electrons more this atom can take
	let electrons_felt = atom.lone_electrons;
	for (let bond of atom.bonds) electrons_felt += bond.order * 2;
	if (atom.atom_object.period == 1) return 2 - electrons_felt;
	if (atom.atom_object.z <= 14) return 8 - electrons_felt;
	return 18 - electrons_felt; //is there any upper limit on the expanded octet rule???
}

function too_many_electrons(atom) {
	let electrons_felt = atom.lone_electrons;
	for (let bond of atom.bonds) electrons_felt += bond.order * 2;
	if (atom.atom_object.period == 1) return electrons_felt > 2;
	if (atom.atom_object.z <= 14) return electrons_felt > 8;
	return electrons_felt > 18; //is there any upper limit on the expanded octet rule???
}

function atom_formal_charge(atom) {
	let immediate_electrons = atom.lone_electrons;
	for (let bond of atom.bonds) immediate_electrons += bond.order;
	return atom.brings - immediate_electrons;
}

//TODO: use this as helper function in others
function atom_electrons_felt(atom) {
	let electrons_felt = atom.lone_electrons;
	for (let bond of atom.bonds) electrons_felt += bond.order * 2;
	return electrons_felt;
}

function lewis_text(molecule = molecule) {
	let text = [];
	for (let atom of molecule.atoms) {
		let str = atom.symbol + " ";
		for (let i = 0; i < atom.lone_electrons - 1; i += 2) str += ":";
		if (atom.lone_electrons % 2) str += ".";
		for (let bond of atom.bonds) {
			str += " " + "*-=≡≣".charAt(bond.order) + bond_other_atom(atom, bond).symbol;
		}
		if (atom_formal_charge(atom) != 0) str += " (charge "+atom_formal_charge(atom)+")";
		text.push(str);
	}
	return text.join("\n");
}

function bond_other_atom(my_atom, bond) {
	if (bond.atoms[0] == my_atom) return bond.atoms[1];
	return bond.atoms[0];
}
function remove_value(arr, val) {return arr.filter(function(e) {return e != val;});}

function are_bonded(atom_a, atom_b) {
	for (let bond of atom_a.bonds) if (bond.atoms.includes(atom_b)) return true;
	return false;
}
function steric_number(atom) {
	return Math.ceil(atom.lone_electrons / 2) + atom.bonds.length;
}