#include "../../src/Forces/Nanopore.h"
#include "../../src/Forces/CylinderConfinement.h"
#include "../../src/Forces/ForceFactory.h"
#include "../../src/Particles/Molecule.h"
#include "../../src/Boxes/CubicBox.h"
#include "../../src/Utilities/Logger.h"

#include <array>
#include <iostream>
#include <random>
#include <stdexcept>

namespace CG = CylindricalGeometry;
int checks = 0;

void require(bool ok, const std::string &label) {
	++checks;
	if(!ok) throw std::runtime_error(label);
}

void near(double got, double expected, double tolerance, const std::string &label) {
	require(std::isfinite(got) && std::abs(got - expected) <= tolerance * (1. + std::abs(expected)),
		label + ": got " + std::to_string(got) + ", expected " + std::to_string(expected));
}

template<typename T> using Evaluator = CG::Evaluation<T> (*)(const CG::Parameters<T> &, T, T, T);

template<typename T>
void gradient_checks(Evaluator<T> evaluate, bool is_pore) {
	std::mt19937 rng(1907);
	std::uniform_real_distribution<double> random(-7., 7.);
	const T delta = sizeof(T) == sizeof(float) ? T(0.002) : T(1.e-5);
	const double tolerance = sizeof(T) == sizeof(float) ? 0.015 : 2.e-7;
	for(int i = 0; i < 2500; ++i) {
		CG::Parameters<T> p{T(3.5), T(2.3), T(3.1), T(i % 2 ? 0.4 : 0.),
			T(1.2), T(-0.7), T(0.3), T(2.), T(-3.), T(4.)};
		std::array<T, 3> point{{T(random(rng)), T(random(rng)), T(random(rng))}};
		const auto local = CG::local_coordinates(p, point[0], point[1], point[2]);
		// The membrane medial axis has no unique derivative. Test it separately
		// below, and keep finite-difference stencils away from it here.
		if(is_pore && (std::abs(local.z) < 20*delta ||
			std::abs((p.radius - local.rho) - (std::abs(local.z) - p.half_length)) < 20*delta)) continue;
		const auto result = evaluate(p, point[0], point[1], point[2]);
		const T forces[] = {result.fx, result.fy, result.fz};
		for(int axis = 0; axis < 3; ++axis) {
			auto plus = point, minus = point;
			plus[axis] += delta;
			minus[axis] -= delta;
			const T numerical = -(evaluate(p, plus[0], plus[1], plus[2]).energy -
				evaluate(p, minus[0], minus[1], minus[2]).energy) / (2*delta);
			near(forces[axis], numerical, tolerance, "force = -gradient(energy)");
		}
		// Reversing the axis cannot change an unoriented membrane/container.
		p.nx = -p.nx; p.ny = -p.ny; p.nz = -p.nz;
		const auto reversed = evaluate(p, point[0], point[1], point[2]);
		near(reversed.energy, result.energy, tolerance, "axis reversal energy");
		near(reversed.fx, result.fx, tolerance, "axis reversal force x");
		near(reversed.fy, result.fy, tolerance, "axis reversal force y");
		near(reversed.fz, result.fz, tolerance, "axis reversal force z");
	}
}

void geometry_cases() {
	CG::Parameters<double> p{10., 2., 3., 0., 0., 0., 0., 0., 0., 1.};
	auto r = CG::cylinder(p, 3., 0., 4.);
	near(r.energy, 10., 1.e-12, "corner adds side and cap energies");
	near(r.fx, -10., 1.e-12, "corner side force");
	near(r.fz, -10., 1.e-12, "corner cap force");
	r = CG::cylinder(p, 0., 0., -4.);
	near(r.fx, 0., 1.e-12, "axis has no radial singularity");
	near(r.fz, 10., 1.e-12, "lower cap points inward");
	r = CG::cylinder(p, 2., 0., 3.);
	near(r.energy, 0., 1.e-12, "closed cylinder surface");
	near(r.fx, 0., 1.e-12, "side force continuous at contact");
	near(r.fz, 0., 1.e-12, "cap force continuous at contact");
	p.particle_radius = 0.5;
	r = CG::cylinder(p, 2., 0., 3.);
	near(r.energy, 2.5, 1.e-12, "finite radius shrinks both dimensions");
	near(r.fx, -5., 1.e-12, "finite radius side");
	near(r.fz, -5., 1.e-12, "finite radius cap");
	// A spherical particle near the pore rim sees Euclidean rim clearance.
	r = CG::nanopore(p, 1.7, 0., 3.4);
	near(r.energy, 0., 1.e-12, "3-4-5 rim contact");
	r = CG::nanopore(p, 1.85, 0., 3.2);
	near(r.energy, 0.3125, 1.e-12, "rim overlap energy");
	near(r.fx, -1.5, 1.e-12, "rim pushes toward pore");
	near(r.fz, 2., 1.e-12, "rim pushes toward reservoir");
	p.particle_radius = 0.;
	r = CG::nanopore(p, 5., 0., 0.);
	near(r.fz, 0., 1.e-12, "midplane does not pick +z");
	r = CG::nanopore(p, 3., 0., 2.);
	near(r.energy, 5., 1.e-12, "medial tie energy");
	near(r.fx, -5., 1.e-12, "medial tie symmetric radial subgradient");
	near(r.fz, 5., 1.e-12, "medial tie symmetric axial subgradient");
	// Legacy point-centre energy and force away from non-differentiable ties.
	for(int i = 0; i < 200; ++i) {
		const double rho = 0.05*i, z = -4.7 + 0.047*i;
		const double d = std::max(0., std::min(rho - p.radius, p.half_length - std::abs(z)));
		r = CG::nanopore(p, rho, 0., z);
		near(r.energy, 0.5*p.stiff*d*d, 1.e-12, "legacy energy");
		if(std::abs((rho-p.radius) - (p.half_length-std::abs(z))) < 1.e-10 || std::abs(z) < 1.e-10) continue;
		if(d > 0. && rho-p.radius < p.half_length-std::abs(z)) {
			near(r.fx, -p.stiff*d, 1.e-12, "legacy radial force");
			near(r.fz, 0., 1.e-12, "legacy radial branch z");
		}
		else if(d > 0.) {
			near(r.fx, 0., 1.e-12, "legacy face branch x");
			near(r.fz, std::copysign(p.stiff*d, z), 1.e-12, "legacy face force");
		}
	}
	// Analytic rotational/translation covariance, independent of the projection code.
	for(auto evaluate : {CG::nanopore<double>, CG::cylinder<double>}) {
		const auto original = evaluate(p, 2.4, 0.5, 2.8);
		auto rotated = p;
		rotated.cx = 1.; rotated.cy = -2.; rotated.cz = 3.;
		rotated.nx = 0.6; rotated.ny = 0.; rotated.nz = 0.8;
		const auto transformed = evaluate(rotated, 1.+0.8*2.4+0.6*2.8, -2.+0.5, 3.-0.6*2.4+0.8*2.8);
		near(transformed.energy, original.energy, 1.e-12, "rotation energy");
		near(transformed.fx, 0.8*original.fx+0.6*original.fz, 1.e-12, "rotation force x");
		near(transformed.fy, original.fy, 1.e-12, "rotation force y");
		near(transformed.fz, -0.6*original.fx+0.8*original.fz, 1.e-12, "rotation force z");
	}
}

template<typename Force>
void invalid_input(const std::string &type, const std::string &length_key, const std::string &extra) {
	input_file inp;
	inp.show_overwrite_warnings = false;
	inp.init_from_string("type=" + type + "\nparticle=all\nstiff=10\nradius=2\n" +
		length_key + "=6\ndir=0,0,1\n" + extra);
	bool rejected = false;
	try { Force force; force.init(inp); }
	catch(const oxDNAException &) { rejected = true; }
	require(rejected, "invalid input rejected: " + type + " " + extra);
}

template<typename Force>
void input_cases(const std::string &type, const std::string &length_key, CubicBox &box) {
	for(const std::string bad : {"stiff=0", "stiff=-1", "stiff=nan", "stiff=inf", "stiff=1x",
		"radius=0", "radius=-1", "radius=1e999", "radius=garbage", "particle_radius=-1",
		"particle_radius=2", "dir=0,0,0", "dir=1,2", "dir=1,2,3,4", "dir=1,2,3x",
		"dir=1,2,3,", "dir=1,inf,0", "center=0,nan,0", "center=0,0,0junk"}) {
		invalid_input<Force>(type, length_key, bad + "\n");
	}
	invalid_input<Force>(type, length_key, length_key + "=0\n");
	invalid_input<Force>(type, length_key, length_key + "=-1\n");
	invalid_input<Force>(type, length_key, length_key + "=nan\n");
	for(const char *settings : {"use_barostat=1\n", "ensemble=NPT\n"}) {
		input_file simulation;
		simulation.init_from_string(settings);
		CONFIG_INFO->sim_input = &simulation;
		invalid_input<Force>(type, length_key, "");
		CONFIG_INFO->sim_input = nullptr;
	}
	input_file inp;
	inp.init_from_string("type=" + type + "\nparticle=all\nstiff=10\nradius=2\n" + length_key +
		"=6\ndir=0,0,1e30\ncenter=1,2,3\ngroup_name=geometry\nid=test_" + type);
	auto molecule = CONFIG_INFO->molecules()[0];
	molecule->invalidate_shiftable();
	require(molecule->shiftable(), "unconstrained molecule is initially shiftable");
	auto *shifted_particle = CONFIG_INFO->particles()[0];
	shifted_particle->pos = LR_vector(4., 2., 3.5);
	shifted_particle->set_pos_shift(1, 0, 0);
	const auto original_absolute = box.get_abs_pos(shifted_particle);
	ForceFactory::instance()->add_force(inp, CONFIG_INFO->particles(), &box);
	require(!molecule->shiftable(), "force attachment invalidates cached shiftability");
	near(shifted_particle->pos.x, original_absolute.x, 1.e-12, "late force attachment restores absolute coordinates for CUDA");
	require(shifted_particle->_pos_shift[0] == 0, "late force attachment clears box shift");
	auto *force = dynamic_cast<Force *>(CONFIG_INFO->forces.back().get());
	require(force != nullptr, "factory registration " + type);
	near(force->_direction.z, 1., 1.e-12, "large axis safely normalised");
	require(force->get_group_name() == "geometry", "force group retained");
	// The force must use absolute coordinates even when a particle has box shifts.
	auto *particle = CONFIG_INFO->particles()[0];
	particle->pos = LR_vector(4., 2., 3.5);
	particle->set_pos_shift(1, 0, 0);
	particle->set_initial_forces(0, &box);
	particle->set_ext_potential(0, &box);
	auto absolute = box.get_abs_pos(particle);
	const auto direct = force->value(0, absolute);
	near(particle->force.x, direct.x, 1.e-12, "particle absolute force");
	near(particle->ext_potential, force->potential(0, absolute), 1.e-12, "particle absolute energy");
	for(auto part : CONFIG_INFO->particles()) part->ext_forces.clear();
	CONFIG_INFO->forces.clear();
}

int main() {
	try {
		Logger::init();
		Logger::instance()->disable_log("cylindrical force tests");
		BaseParticle p0, p1;
		p0.index = 0; p1.index = 1;
		p0.strand_id = p1.strand_id = 0;
		std::vector<BaseParticle *> particles{&p0, &p1};
		std::vector<std::shared_ptr<Molecule>> molecules{std::make_shared<Molecule>()};
		for(auto particle : particles) molecules[0]->add_particle(particle);
		ConfigInfo::init(&particles, &molecules);
		CubicBox box;
		box.init(30., 30., 30.);
		CONFIG_INFO->box = &box;
		geometry_cases();
		gradient_checks<double>(CG::nanopore<double>, true);
		gradient_checks<double>(CG::cylinder<double>, false);
		gradient_checks<float>(CG::nanopore<float>, true);
		gradient_checks<float>(CG::cylinder<float>, false);
		input_cases<Nanopore>("nanopore", "thickness", box);
		input_cases<CylinderConfinement>("cylinder_confinement", "height", box);
		invalid_input<CylinderConfinement>("cylinder_confinement", "height", "height=1\nparticle_radius=0.5\n");
		ConfigInfo::clear();
		std::cout << "Passed " << checks << " cylindrical force checks (geometry, gradients, float/double, input, factory, absolute positions).\n";
		return 0;
	}
	catch(const std::exception &error) {
		std::cerr << "FAILED after " << checks << " checks: " << error.what() << '\n';
		return 1;
	}
}
