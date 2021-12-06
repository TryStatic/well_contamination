#include <iostream>
#include <vector>
#include <cmath>

#pragma warning(disable:4996)

// Σταθερές
#define PI				3.14159265358979323846
#define Q				0.03
#define K				0.0001
#define A				40
#define N				0.2
#define VE				0.000005
#define HALT_DISTANCE	50.0
#define BACKTRACK_STEPS 9

// Δομή για τα σημεία
struct point_2d
{
	double x;
	double y;

	point_2d(const double x, const double y)
	{
		this->x = x;
		this->y = y;
	}
};

// Δομή για τα σωματίδια
struct particle
{
	point_2d position{0,0};
	int halted_step;

	particle()
	{
		position = point_2d(0, 0);
		halted_step = -1;
	}

	particle(const double x, const double y)
	{
		position = point_2d(x, y);
		halted_step = -1;
	}

	particle(const point_2d& point)
	{
		position = point_2d(point.x, point.y);
		halted_step = -1;
	}
};

// Δομη ορισμού ενός πηγαδιού
struct well
{
	point_2d center;
	double radius;
	
	well(const point_2d& center, const double radius)	
		: center(center), radius(radius)
	{
	}

	// Μέθοδος υπολογισμού σημείων περιμετρικά του πηγαδιού
	std::vector<point_2d> generate_particles(const int particles_count) const
	{
		std::vector<point_2d> points;

		for (int i = 0; i != particles_count; ++i)
		{
			const double angle_rad = 2 * i * PI / particles_count;

			double x = this->radius * cos(angle_rad) + this->center.x;
			double y = this->radius * sin(angle_rad) + this->center.y;

			points.emplace_back(x, y);
		}

		return points;
	}
};

// Συνάρτηση υπολογισμού απόστασης δύο σημείων
double get_distance(const point_2d& point1, const point_2d& point2)
{
	return sqrt(pow(point1.x - point2.x, 2) + pow(point1.y - point2.y, 2));
}

// Συνάρτηση υπολογισμού σωματιδίων μόλυνσης
void track_contamination( // Παράμετροι:
	FILE* f, // Το αρχείο για το output
	const well well1, // To πρώτο πηγάδι
	const well well2, // To δεύτερο πηγάδι
	const std::vector<point_2d>& particle_positions, // Λίστα με τα αρχικά σημεία
	const int steps, // Αρ. βημάτων
	const int deltaT, // το dt
	const bool inverse = false // Εάν υπολογίζουμε ανάποδα ή οχι
	)
{
	// Λίστα για τα σωματίδια
	std::vector<particle> particles;

	// Δέσμευση χώρου (μνήμης) για τα σωματίδια
	particles.reserve(particle_positions.size());

	// Για κάθε σημείο
	for (point_2d p : particle_positions)
	{
		// Αρχικοποιούμε ενα particle / σωματήδιο
		particles.emplace_back(p);
	}

	std::vector<std::vector<particle>> particles_memory;
	particles_memory.emplace_back(particles);

	// Εκτύπωση πληροφοριών για τήν τρέχουσα προσομοίωση
	fprintf(f, "Tracking contamination for %3d particles and %3d time steps. Reverse tracking: %s\n", particles.size(), steps, inverse ? "Enabled" : "Disabled");
	fprintf(f, "Well 1 | Center = (%10.4f, %10.4f), Radius=%10.4f\n", well1.center.x, well1.center.y, well1.radius);
	fprintf(f, "Well 2 | Center = (%10.4f, %10.4f), Radius=%10.4f\n", well2.center.x, well2.center.y, well2.radius);

	// Εκτύπωση των αρχικών σημείων
	fprintf(f, "Initial particle positions: \n");
	for(int i = 0; i != particles.size(); ++i)
	{
		fprintf(f, "particle=%3d | x=%10.4f | Y=%10.4f\n", i, particles[i].position.x, particles[i].position.y);
	}

	// Για κάθε βημα
	for (int i = 0; i != steps; i++)
	{
		fprintf(f, "STEP %3d | DAY %4d\n", i, (i+1) * deltaT / (60 * 60 * 24));
		// + Για κάθε σωματίδιο στο τρέχων βήμα
		for (unsigned int j = 0; j != particles.size(); j++)
		{
			// Υπολογισμός ταχυτήτων
			
			const double c = Q / (2 * A * N * PI);

			const double particle_x = particles[j].position.x;
			const double particle_y = particles[j].position.y;

			const double well1_x = well1.center.x;
			const double well1_y = well1.center.y;

			const double well2_x = well2.center.x;
			const double well2_y = well2.center.y;

			const double xxa = particle_x - well1_x;
			const double yya = particle_y - well1_y;

			const double yyb = particle_y - well2_y;
			const double xxb = particle_x - well2_x;

			const double p1_x = (xxa) / (pow(xxa, 2) + pow(yya, 2));
			const double p2_x = (xxb) / (pow(xxb, 2) + pow(yyb, 2));


			const double p1_y = (yya) / (pow(xxa, 2) + pow(yya, 2));
			const double p2_y = (yyb) / (pow(xxb, 2) + pow(yyb, 2));

			double vx = c * (p1_x + p2_x);
			double vy = c * (p1_y + p2_y);

			double inverse_modifier = 1.0;

			// Αντιστρέφουμε πρόσημα εαν κανουμε ανάποδη προσομοιωση
			if(inverse)
			{
				vx *= -1;
				vy *= -1;
				inverse_modifier = -1;
			}

			const double new_position_x = particle_x + (vx * deltaT);
			const double new_position_y = particle_y + (vy * deltaT) + (VE * deltaT * inverse_modifier);

			// Εάν σε κάποιο προιγούμενο βήμα έχουμε υπολογίσει οτι το συγκεκριμένο σωματίδιο έχει σταματήσει
			// Εκτυπώνουμε απλά σε ποιο βήμα σταματήσαμε.
			if(particles[j].halted_step != -1)
			{
				fprintf(f, "particle=%3d | newX=%10.4f | newY=%10.4f | HALT AT STEP=%2d\n", j, particles[j].position.x, particles[j].position.y, particles[j].halted_step);
				continue;
			}

			// Μεταβλητή ελέγχου (εαν κάποιο σωματήδιο πάει κοντά σε κάποιο απο τα πηγάδια θα πάρει τον αριθμό του πηγαδιού... 1 ή 2)
			int getting_close_to_well = -1;

			// Check-aroume αποσταση particle + πηγάδι 1
			if (get_distance(well1.center, point_2d(new_position_x, new_position_y)) < HALT_DISTANCE)
			{
				fprintf(f, "Particle=%3d is halted at position (%10.4f, %10.4f) due to being close to well 1\n", j, particle_x, particle_y);
				particles[j].halted_step = i;
				getting_close_to_well = 1;
			}

			// Check-aroume αποσταση particle + πηγάδι 2
			if (get_distance(well2.center, point_2d(new_position_x, new_position_y)) < HALT_DISTANCE)
			{
				fprintf(f, "Particle=%3d is halted at position (%10.4f, %10.4f) due to being close to well 2\n", j, particle_x, particle_y);
				particles[j].halted_step = i;
 				getting_close_to_well = 2;
			}

			// Εαν η μεταβλητή μας είναι διαφορετική από το -1, συμαίνει οτι το συγκεκριμένο σωματήδιο (j) στο συγκεκριμένο βήμα (i)
			// κοντεύει να μπεί στο πηγάδι (dist < 50m). Σταματάμε εδώ και κάνουμε backtrack στην μνήμη να βρούμε την τοποθεσία του
			// 6 μήνες πριν.
			if(getting_close_to_well != -1)
			{
				// Πρέπει να βρούμε ποιο βήμα αντιστοιχεί σε 6 μήνες πρίν.
				// Αν 1 βήμα είναι 20 ήμερες τότε 9 βήματα θα έιναι 180 ήμερες ή αλλιώς 6 μήνες.
				// Αρά θέλουμε να βρούμε που ήταν το σωματήδιο 9 βήματα πριν (i - 9)

				int backtrack_step = i - 9;

				// Sanity check
				if (backtrack_step < 0) backtrack_step = 0;

				double pos_x = particles_memory[backtrack_step][j].position.x;
				double pos_y = particles_memory[backtrack_step][j].position.y;
				
				fprintf(f, "ALERT for particle %2d at position (%10.4f, %10.4f) would contamine well %2d within 6 months.\n", j, pos_x, pos_y, getting_close_to_well);
				return;
			}

			// Ενημερώνουμε συνεταγμένες στον πίνακα μας / λίστα / vector
			particles[j].position.x = particle_x + (vx * deltaT);
			particles[j].position.y = particle_y + (vy * deltaT) + (VE * deltaT * inverse_modifier);

			// Βάζουμε τις συντεταγμένες στην μνήμη.
			particles_memory.emplace_back(particles);
			
			// Εκτυπώνουμε
			fprintf(f, "particle=%3d | newX=%10.4f | newY=%10.4f\n", j, particles[j].position.x, particles[j].position.y);
		}
	}
}

int main()
{
	FILE* file = fopen("output.txt", "w");

	if(!file)
	{
		std::cout << "Could not open output.txt file for writing." << std::endl;
		return 0;
	}
	
	const int steps = 92;
	const int particle_count = 10;

	// Wells
	const well well1 = well(point_2d(1000, 150), 50);
	const well well2 = well(point_2d(300, 100), 50);

	const std::vector<point_2d> particles1 =
	{
		point_2d(1046.19, 169.13),
		point_2d(1035.36, 185.36),
		point_2d(1019.13, 196.19),
		point_2d(1000.0, 200.0),
		point_2d(980.87, 196.19),
		point_2d(964.64, 185.36),
		point_2d(953.81, 169.13),
		point_2d(950.0, 150.0),
		point_2d(953.81, 130.87),
		point_2d(964.64, 114.64),
		point_2d(980.87, 103.81),
		point_2d(1000.0, 100.0),
		point_2d(1019.13, 103.81),
		point_2d(1035.36, 114.64),
		point_2d(1046.19, 130.87),
		point_2d(1050.0, 150.0),
	};

	track_contamination(file, well1, well2, particles1, 92, 1728000);


	const std::vector<point_2d> particles2 =
	{
		point_2d(346.19,	119.13),
		point_2d(335.36,	135.36),
		point_2d(319.13,	146.19),
		point_2d(300.00,	150.00),
		point_2d(280.87,	146.19),
		point_2d(264.64,	135.36),
		point_2d(253.81,	119.13),
		point_2d(250.00,	100.00),
		point_2d(253.81,	80.87),
		point_2d(264.64,	64.64),
		point_2d(280.87,	53.81),
		point_2d(300.00,	50.00),
		point_2d(319.13,	53.81),
		point_2d(335.36,	64.64),
		point_2d(346.19,	80.87),
		point_2d(350.00,	100.00),
	};

	track_contamination(file, well1, well2, particles2, 92, 1728000);

	
	const std::vector<point_2d> particles3 =
	{
		point_2d(80.00, 1195.00),
		point_2d(85.00, 1195.00),
		point_2d(90.00, 1195.00),
		point_2d(95.00, 1195.00),
		point_2d(100.00, 1195.00),
		point_2d(105.00, 1195.00),
		point_2d(110.00, 1195.00),
		point_2d(115.00, 1195.00),
		point_2d(120.00, 1195.00),
		point_2d(125.00, 1195.00),
	};

	track_contamination(file, well1, well2, particles3, 92, 1728000, true);


	const std::vector<point_2d> particles4 =
	{
		point_2d(280.00,	695.00),
		point_2d(285.00,	695.00),
		point_2d(290.00,	695.00),
		point_2d(295.00,	695.00),
		point_2d(300.00,	695.00),
		point_2d(305.00,	695.00),
		point_2d(310.00,	695.00),
		point_2d(315.00,	695.00),
		point_2d(320.00,	695.00),
		point_2d(325.00,	695.00)
	};

	track_contamination(file, well1, well2, particles4, 92, 1296000, true);

	const std::vector<point_2d> particles5 =
	{
		point_2d(780,	1095),
		point_2d(785,	1095),
		point_2d(790,	1095),
		point_2d(795,	1095),
		point_2d(800,	1095),
		point_2d(805,	1095),
		point_2d(810,	1095),
		point_2d(815,	1095),
		point_2d(820,	1095),
		point_2d(825,	1095),
	};

	track_contamination(file, well1, well2, particles5, 92, 1728000, true);

	const std::vector<point_2d> particles6 =
	{
		point_2d(1180.00, 995.00),
		point_2d(1185.00, 995.00),
		point_2d(1190.00, 995.00),
		point_2d(1195.00, 995.00),
		point_2d(1200.00, 995.00),
		point_2d(805.00, 995.00),
		point_2d(1210.00, 995.00),
		point_2d(1215.00, 995.00),
		point_2d(1220.00, 995.00),
		point_2d(1225.00, 995.00),
	};

	track_contamination(file, well1, well2, particles6, 92, 1728000, true);

	fclose(file);
}