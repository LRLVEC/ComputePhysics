#include <_BLAS.h>
#include <_Time.h>

using namespace BLAS;

struct Thomson
{
	static constexpr double answers[66]
	{
		0,
		0,
		0.500000000,
		1.732050808,
		3.674234614,
		6.474691495,
		9.985281374,
		14.452977414,
		19.675287861,
		25.759986531,
		32.716949460,
		40.596450508,
		49.165253058,
		58.853230612,
		69.306363297,
		80.670244114,
		92.911655303,
		106.050404829,
		120.084467447,
		135.089467557,
		150.881568334,
		167.641622399,
		185.287536149,
		203.930190663,
		223.347074052,
		243.812760299,
		265.133326317,
		287.302615033,
		310.491542358,
		334.634439920,
		359.603945904,
		385.530838063,
		412.261274651,
		440.204057448,
		468.904853281,
		498.569872491,
		529.122408375,
		560.618887731,
		593.038503566,
		626.389009017,
		660.675278835,
		695.916744342,
		732.078107544,
		769.190846459,
		807.174263085,
		846.188401061,
		886.167113639,
		927.059270680,
		968.713455344,
		1011.557182654,
		1055.182314726,
		1099.819290319,
		1145.418964319,
		1191.922290416,
		1239.361474729,
		1287.772720783,
		1337.094945276,
		1387.383229253,
		1438.618250640,
		1490.773335279,
		1543.830400976,
		1597.941830199,
		1652.909409898,
		1708.879681503,
		1765.802577927
	};
	unsigned long long num;
	vec pos;
	vec g0;
	double answer;

	Thomson(unsigned long long _num)
		:
		num(_num),
		pos(_num * 2, false),
		g0(_num * 2, false)
	{
	}
	double rij(double theta0, double phi0, double theta1, double phi1)
	{
		double st12(sin((theta0 - theta1) / 2));
		double sp12(sin((phi0 - phi1) / 2));
		return 2 * sqrt(st12 * st12 + sin(theta0) * sin(theta1) * sp12 * sp12);
	}
	double psi(unsigned long long c0, double theta0, double phi0)
	{
		double p(0);
		for (unsigned long long c1(0); c1 < c0; ++c1)
			p += 1.0 / rij(theta0, phi0, pos[2 * c1], pos[2 * c1 + 1]);
		for (unsigned long long c1(c0 + 1); c1 < num; ++c1)
			p += 1.0 / rij(theta0, phi0, pos[2 * c1], pos[2 * c1 + 1]);
		return p;
	}
	double psi(vec const& p)
	{
		double ss(0);
		for (unsigned long long c0(1); c0 < num; ++c0)
		{
			double theta0(p.data[2 * c0]);
			double phi0(p.data[2 * c0 + 1]);
			for (unsigned long long c1(0); c1 < c0; ++c1)
				ss += 1 / rij(theta0, phi0, pos[2 * c1], pos[2 * c1 + 1]);
		}
		return ss;
	}
	void initPos(std::uniform_real_distribution<double>& rd, std::mt19937& mt)
	{
		for (unsigned long long c0(0); c0 < num; ++c0)
		{
			pos[2 * c0] = acos(2 * rd(mt) - 1);
			pos[2 * c0 + 1] = 2 * Pi * rd(mt);
		}
	}
	void gradient(vec const& p, vec& g)
	{
		constexpr double h = 1e-5;
		constexpr double rh = 1 / (2.0 * h);
		g = 0;
		//vec gg(2 * num, true);
		for (unsigned long long c0(0); c0 < num; ++c0)
		{
			double theta0 = p.data[2 * c0];
			double phi0 = p.data[2 * c0 + 1];
			double gt0, gt1, gp(0);
			//g[2 * c0] = rh * (psi(c0, theta0 + h, phi0) - psi(c0, theta0 - h, phi0));
			//g[2 * c0 + 1] = rh * (psi(c0, theta0, phi0 + h) - psi(c0, theta0, phi0 - h));
			for (unsigned long long c1(0); c1 < c0; ++c1)
			{
				double theta1 = p.data[2 * c1];
				double phi1 = p.data[2 * c1 + 1];
				double s0(sin(theta0)), s1(sin(theta1));
				double st12(sin((theta0 - theta1) / 2));
				double sp12(sin((phi0 - phi1) / 2));
				sp12 *= sp12;
				double ss0(s0 * sp12);
				double ss1(s1 * sp12);
				double rr(2 * sqrt(st12 * st12 + s0 * ss1));
				double r3(rr * rr * rr);
				double s01(sin(theta0 - theta1));
				gt0 = (s01 + 2 * cos(theta0) * ss1) / r3;
				gt1 = (-s01 + 2 * cos(theta1) * ss0) / r3;
				gp = (s0 * s1 * sin(phi0 - phi1)) / r3;
				g[2 * c0] -= gt0;
				g[2 * c1] -= gt1;
				g[2 * c0 + 1] -= gp;
				g[2 * c1 + 1] += gp;
			}
		}
		//gg -= g;
		//::printf("gg: %.10e\n", gg.norm2());
	}
	void naive()
	{
		::printf("%.5e\n", psi(pos) - 1765.802577927);
		for (unsigned long long c0(0); c0 < 2000; ++c0)
		{
			gradient(pos, g0);
			pos.fmadd(-2e-3, g0);
		}
		::printf("%.5e\n", psi(pos) - 1765.802577927);
	}
	double conjugateGradientDavidon()
	{
		vec pos1(num * 2, false);
		vec postp[2]{ vec(num * 2, false),  vec(num * 2, false) };
		vec g1(num * 2, false);
		vec gtp[2]{ vec(num * 2, false), vec(num * 2, false) };
		vec y(num * 2, false);
		vec d(num * 2, false);
		vec answer(100, true);
		double tau((sqrt(5) - 1) / 2);
		double tau1(1 - tau);
		double alpha;
		double beta;
		double psi0(psi(pos));
		double psi1;
		gradient(pos, g0);
		double gNorm(g0.norm2Square()), gNorm1;
		double gm(gNorm);
		double dg0(-gNorm), dg1;
		unsigned long long gidx(0);//0: g0, 1: g1, 2: gtp[0], 3: gtp[1]
		unsigned long long cnt(0);
		unsigned long long an(0);
		while (gNorm > 1e-18)
		{
			beta = sqrt(gNorm) * 1e-6;
			if (dg0 > 0 || cnt % (2 * num) == 0)
			{
				d = g0;
				d.neg();
				dg0 = -gNorm;
			}
			for (;;)
			{
				pos1.fmadd(beta, d, pos);
				gradient(pos1, g1);
				dg1 = (d, g1);
				psi1 = psi(pos1);
				if (dg1 > 0)break;
				beta *= 2;
			}
			while (-dg0 > 0.5 * gNorm)
			{
				double e(3 * (psi0 - psi1) / beta + dg0 + dg1);
				double h(sqrt(e * e - dg0 * dg1));
				alpha = beta;
				double kk(dg0 + dg1 + 2 * e);
				beta *= (dg0 + e + h) / kk;
				pos1.fmadd(beta, d, pos);
				gradient(pos1, g1);
				double dg = (d, g1);
				if ((alpha - beta) / alpha < 0.01)
				{
					double ds(2 * (alpha - beta));
					double dgtp0(dg1), dgtp1;
					double at(0);
					bool id(false);
					for (;;)
					{
						if (ds > alpha)break;
						postp[id].fmadd(alpha - ds, d, pos);
						gradient(postp[id], gtp[id]);
						dgtp1 = dgtp0;
						dgtp0 = (d, gtp[id]);
						if (dgtp0 < 0)break;
						at = ds;
						ds *= 2;
						id = !id;
					}
					if (ds < alpha)
					{
						pos = postp[id];
						psi0 = psi(pos);
						dg0 = dgtp0;
						gidx = 2 + unsigned long long(id);
					}
					if (at != 0)
					{
						pos1 = postp[!id];
						psi1 = psi(pos1);
						dg1 = dgtp1;
					}
					beta = ds - at;
					continue;
				}
				else if (beta / alpha < 0.01)
				{
					double ds(beta * 2);
					double dgtp0, dgtp1(dg0);
					double at(0);
					bool id(false);
					for (;;)
					{
						if (ds > alpha)break;
						postp[id].fmadd(ds, d, pos);
						gradient(postp[id], gtp[id]);
						dgtp0 = dgtp1;
						dgtp1 = (d, gtp[id]);
						if (dgtp1 > 0)break;
						at = ds;
						ds *= 2;
						id = !id;
					}
					if (ds < alpha)
					{
						pos1 = postp[id];
						dg1 = dgtp1;
						psi1 = psi(pos1);
					}
					if (at != 0)
					{
						pos = postp[!id];
						psi0 = psi(pos);
						dg0 = dgtp0;
						gidx = 2 + unsigned long long(!id);
					}
					beta = ds - at;
					continue;
				}
				if (dg == 0)break;
				else if (dg > 0) { psi1 = psi(pos1); dg1 = dg; }
				else { psi0 = psi(pos1); dg0 = dg; pos = pos1; beta = alpha - beta; gidx = 1; }
			}
			switch (gidx)
			{
			case 2: g1 = gtp[0]; break;
			case 3: g1 = gtp[1]; break;
			}
			y = g1;
			y -= g0;
			g0 = g1;
			gNorm1 = gNorm;
			gNorm = g0.norm2Square();
			double dy((d, y));
			double y2(y.norm2Square());
			y.fmadd(-2.0 * y2 / dy, d);
			double dBeta((y, g0) / dy);
			d *= dBeta;
			d -= g0;
			dg0 = (d, g0);
			//::printf("%.10e\t%.10e\t%.10e\n", psi0 - answers[num], dg0, gNorm);
			if (gNorm < 1e-5 && abs(psi0 - answers[num])>2e-4)
			{
				answer[cnt % answer.dim] = psi0;
				++an;
				break;
			}
			if (gNorm < (num > 20 ? 1e-13 : 1e-16))
			{
				answer[cnt % answer.dim] = psi0;
				++an;
			}
			++cnt;
		}
		double sm(answer.norm1() / (an < answer.dim ? an : answer.dim));
		//double delta(sm - answers[num]);
		//::printf("%llu: %.3e\n", an, delta);
		return sm;
	}
	double run(std::uniform_real_distribution<double>& rd, std::mt19937& mt)
	{
		do
		{
			initPos(rd, mt);
			answer = conjugateGradientDavidon();
		} while (abs(answer - answers[num]) > 1e-8);
		return answer;
	}
};

struct Icosahedron
{
	Thomson tms;

	Icosahedron(std::uniform_real_distribution<double>& rd, std::mt19937& mt)
		:
		tms(12)
	{
		tms.run(rd, mt);
	}
	double check()
	{
		vec rij(tms.num * (tms.num - 1) / 2, false);
		double min(2);
		for (unsigned long long c0(1); c0 < tms.num; ++c0)
		{
			double theta0(tms.pos[2 * c0]), phi0(tms.pos[2 * c0 + 1]);
			for (unsigned long long c1(0); c1 < c0; ++c1)
			{
				unsigned long long idx(((c0 - 1) * c0) / 2 + c1);
				rij[idx] = tms.rij(theta0, phi0, tms.pos[2 * c1], tms.pos[2 * c1 + 1]);
				if (min > rij[idx])min = rij[idx];
			}
		}
		min += 1e-5;
		unsigned long long minNum(0);
		for (unsigned long long c0(0); c0 < rij.dim; ++c0)
			if (rij[c0] < min)minNum++;
		vec mins(minNum, false);
		minNum = 0;
		for (unsigned long long c0(0); c0 < rij.dim; ++c0)
			if (rij[c0] < min)mins[minNum++] = rij[c0];
		double minAverage(mins.norm1() / mins.dim);
		mins -= minAverage;
		return mins.norm2();
	}
	void ahh()
	{

	}
};

int main()
{
	std::mt19937 mt(time(nullptr));
	std::uniform_real_distribution<double> rd(0, 1.0);
	std::uniform_int_distribution<unsigned long long> rduint(1, 10);
	Timer timer;

	/*timer.begin();
	for (unsigned long long c0(2); c0 < 65; ++c0)
	{
		Thomson tms(c0);
		tms.run(rd, mt);
		::printf("%llu:\t%.10f\t%.3e\n", c0, tms.answer, tms.answer - tms.answers[c0]);
	}
	timer.end();
	timer.print();*/

	Icosahedron ico(rd, mt);
	::printf("%.5e\n", ico.check());
}