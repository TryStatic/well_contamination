#include <iostream>
#include <vector>
#include <cmath>

#define PI 3.14159265358979323846
#define Q 0.03
#define K 0.0001
#define A 40
#define N 0.2
#define VE 0.000005

struct point_2d
{
	double x;
	double y;

	point_2d(const double x, const double y)
	{
		this->x = x;
		this->y = y;
	}

	void print() const
	{
		std::cout << this->x << " " << this->y << std::endl;
	}
};

class well
{
public:
	well(const point_2d& center, const double radius)	
		: center_(center),
		  radius_(radius)
	{
	}

	std::vector<point_2d> generate_particles(const int particles_count) const
	{
		std::vector<point_2d> points;

		for (int i = 0; i != particles_count; ++i)
		{
			const double angle_rad = 2 * i * PI / particles_count;

			double x = this->radius_ * cos(angle_rad) + this->center_.x;
			double y = this->radius_ * sin(angle_rad) + this->center_.y;

			points.emplace_back(x, y);
		}

		return points;
	}

	point_2d center_;
	double radius_;
};

void print_vector(std::vector<point_2d>& vec)
{
	for (point_2d& p : vec)
	{
		std::cout << p.x << " " << p.y << std::endl;
	}
	std::cout << std::endl;
}
void track_contamination(
	const well well1, 
	const well well2, 
	std::vector<point_2d> particles, 
	const int steps, 
	const int deltaT, 
	const bool inverse = false)
{
	// Foreach time step
	std::cout << "Tracking contamination for " << particles.size() << " particles and " << steps << " time steps." << std::endl;
	std::cout << "Well 1: " << "(" << well1.center_.x << ", " << well1.center_.y << "), radius=" << well1.radius_ << std::endl;
	std::cout << "Well 2: " << "(" << well2.center_.x << ", " << well2.center_.y << "), radius=" << well2.radius_ << std::endl << std::endl;

	std::cout << "Initial particle positions: " << std::endl;
	print_vector(particles);
	
	for (int i = 0; i != steps; i++)
	{
		std::cout << "STEP " << i+1 << std::endl;
		// For each particle
		for (int j = 0; j != particles.size(); j++)
		{
			const double c = Q / (2 * A * N * PI);

			const double particle_x = particles[j].x;
			const double particle_y = particles[j].y;

			const double well1_x = well1.center_.x;
			const double well1_y = well1.center_.y;

			const double well2_x = well2.center_.x;
			const double well2_y = well2.center_.y;

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

			if(inverse)
			{
				vx *= -1;
				vy *= -1;
			}
			
			particles[j].x = particle_x + (vx * deltaT);
			particles[j].y = particle_y + (vy * deltaT) + (VE * deltaT);

			std::cout << "oldX=" << particle_x << "\toldY=" << particle_y << "\tVX=" << vx << "\tVY=" << vy << "\tnewX=" << particles[j].x << "\tnewY=" << particles[j].y << std::endl;
		}
	}
}

int main()
{
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

	track_contamination(well1, well2, particles1, 92, 1728000);


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

	track_contamination(well1, well2, particles2, 72, 1296000);

	
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

	track_contamination(well1, well2, particles3, 92, 1728000, true);


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

	track_contamination(well1, well2, particles4, 92, 1296000, true);

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

	track_contamination(well1, well2, particles5, 92, 1728000, true);

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

	track_contamination(well1, well2, particles6, 92, 1728000, true);
}