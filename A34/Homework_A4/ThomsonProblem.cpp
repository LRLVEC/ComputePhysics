#include <_BLAS.h>
#include <_Time.h>

using namespace BLAS;

struct Thomson
{
	static constexpr double answers[101]
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
		1765.802577927,
		1823.667960264,
		1882.441525304,
		1942.122700406,
		2002.874701749,
		2064.533483235,
		2127.100901551,
		2190.649906425,
		2255.001190975,
		2320.633883745,
		2387.072981838,
		2454.369689040,
		2522.674871841,
		2591.850152354,
		2662.046474566,
		2733.248357479,
		2805.355875981,
		2878.522829664,
		2952.569675286,
		3027.528488921,
		3103.465124431,
		3180.361442939,
		3258.211605713,
		3337.000750014,
		3416.720196758,
		3497.439018625,
		3579.091222723,
		3661.713699320,
		3745.291636241,
		3829.844338421,
		3915.309269620,
		4001.771675565,
		4089.154010060,
		4177.533599622,
		4266.822464156,
		4357.139163132,
		4448.350634331,
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
	double conjugateGradientDavidon(double eps)//eps is for gNorm^2
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
		bool succeed(true);
		while (gNorm > eps)
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
				succeed = false;
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
		//if (succeed)::printf("%llu\t", cnt);
		return sm;
	}
	double run(double eps, std::uniform_real_distribution<double>& rd, std::mt19937& mt)
	{
		do
		{
			initPos(rd, mt);
			answer = conjugateGradientDavidon(eps);
		} while (abs(answer - answers[num]) > 1e-8);
		return answer;
	}
};

struct Vibration
{
	Thomson tms;
	mat U;
	vec freqs;
	vec fn;
	vec fs2r;
	unsigned long long fs;

	Vibration(unsigned long long _num, std::uniform_real_distribution<double>& rd, std::mt19937& mt)
		:
		tms(_num),
		U(2 * _num, 2 * _num, false),
		freqs(2 * _num, false),
		fn(2 * _num, true),
		fs2r(2 * _num, true),
		fs(0)
	{
		tms.run(1e-26, rd, mt);
		getUmat(tms.pos, U);
	}
	void getUmat(vec const& p, mat& u)//check this...
	{
		//make sure that mat is in the right format: p.dim square mat
		u.clear();
		for (unsigned long long cy0(1); cy0 < tms.num; ++cy0)
		{
			double theta0 = p.data[2 * cy0];
			double phi0 = p.data[2 * cy0 + 1];
			double dg;
			for (unsigned long long cy1(0); cy1 < cy0; ++cy1)
			{
				double theta1 = p.data[2 * cy1];
				double phi1 = p.data[2 * cy1 + 1];
				double t01(theta0 - theta1);
				double p01(phi0 - phi1);
				double s0(sin(theta0)), s1(sin(theta1));
				double c0(cos(theta0)), c1(cos(theta1));
				double s0s1(s0 * s1);
				double st01(sin(t01)), ct01(cos(t01));
				double sp01(sin(p01)), cp01(cos(p01));
				double cp01m(1 - cp01);
				double pt0(s1 * c0 * cp01m + st01);
				double pt1(s0 * c1 * cp01m - st01);
				double pp0(s0s1 * sp01);//pp1 = -pp0
				double ptt00(ct01 - s0s1 * cp01m);//ptt11 = ptt00
				double ptt01(c0 * c1 * cp01m - ct01);//ptt01 = ptt10
				double ppp00(cp01 * s0s1);//ppp01 = -ppp00, ppp10 = -ppp00, ppp11 = ppp00
				double ptp00(c0 * s1 * sp01);//ptp01 = -ptp00
				double ptp10(c1 * s0 * sp01);//ptp11 = -ptp10
				double r(2 * (1 - ptt00));
				double rn52(pow(r, 2.5));

				dg = (3 * pt0 * pt0 - r * ptt00) / rn52;
				u(2 * cy0, 2 * cy0) += dg;

				dg = (3 * pt1 * pt1 - r * ptt00) / rn52;
				u(2 * cy1, 2 * cy1) += dg;

				dg = (3 * pt0 * pt1 - r * ptt01) / rn52;
				u(2 * cy0, 2 * cy1) += dg;
				u(2 * cy1, 2 * cy0) += dg;

				dg = (3 * pt0 * pp0 - r * ptp00) / rn52;
				u(2 * cy0, 2 * cy0 + 1) += dg;
				u(2 * cy0, 2 * cy1 + 1) -= dg;
				u(2 * cy0 + 1, 2 * cy0) += dg;
				u(2 * cy1 + 1, 2 * cy0) -= dg;

				dg = (3 * pt1 * pp0 - r * ptp10) / rn52;
				u(2 * cy1, 2 * cy0 + 1) += dg;
				u(2 * cy1, 2 * cy1 + 1) -= dg;
				u(2 * cy0 + 1, 2 * cy1) += dg;
				u(2 * cy1 + 1, 2 * cy1) -= dg;

				dg = (3 * pp0 * pp0 - r * ppp00) / rn52;
				u(2 * cy0 + 1, 2 * cy0 + 1) += dg;
				u(2 * cy0 + 1, 2 * cy1 + 1) -= dg;
				u(2 * cy1 + 1, 2 * cy0 + 1) -= dg;
				u(2 * cy1 + 1, 2 * cy1 + 1) += dg;
			}
		}
		/*for (unsigned long long cy0(1); cy0 < tms.num; ++cy0)
		{
			constexpr double h = 1e-5;
			constexpr double rh = 1 / (h * h);
			constexpr double r4h = 1 / (4 * h * h);
			double theta0 = p.data[2 * cy0];
			double phi0 = p.data[2 * cy0 + 1];
			double theta0p(theta0 + h);
			double theta0n(theta0 - h);
			double phi0p(phi0 - h);
			double phi0n(phi0 - h);
			double dg;
			for (unsigned long long cy1(0); cy1 < cy0; ++cy1)
			{
				double theta1 = p.data[2 * cy1];
				double phi1 = p.data[2 * cy1 + 1];
				double theta1p(theta1 + h);
				double theta1n(theta1 - h);
				double phi1p(phi1 - h);
				double phi1n(phi1 - h);
				double r0(1 / tms.rij(theta0, phi0, theta1, phi1));
				double rr[4];
				rr[0] = 1 / tms.rij(theta0p, phi0, theta1, phi1);
				rr[1] = 1 / tms.rij(theta0n, phi0, theta1, phi1);
				dg = rh * (rr[0] + rr[1] - 2 * r0);
				tst(2 * cy0, 2 * cy0) += dg;
				tst(2 * cy1, 2 * cy1) += dg;

				rr[0] = 1 / tms.rij(theta0p, phi0, theta1p, phi1);
				rr[1] = 1 / tms.rij(theta0n, phi0, theta1p, phi1);
				rr[2] = 1 / tms.rij(theta0p, phi0, theta1n, phi1);
				rr[3] = 1 / tms.rij(theta0n, phi0, theta1n, phi1);
				dg = r4h * (rr[0] - rr[1] - rr[2] + rr[3]);
				tst(2 * cy0, 2 * cy1) += dg;
				tst(2 * cy1, 2 * cy0) += dg;

				rr[0] = 1 / tms.rij(theta0p, phi0p, theta1, phi1);
				rr[1] = 1 / tms.rij(theta0n, phi0p, theta1, phi1);
				rr[2] = 1 / tms.rij(theta0p, phi0n, theta1, phi1);
				rr[3] = 1 / tms.rij(theta0n, phi0n, theta1, phi1);
				dg = r4h * (rr[0] - rr[1] - rr[2] + rr[3]);
				tst(2 * cy0, 2 * cy0 + 1) += dg;
				tst(2 * cy0, 2 * cy1 + 1) -= dg;
				tst(2 * cy0 + 1, 2 * cy0) += dg;
				tst(2 * cy1 + 1, 2 * cy0) -= dg;

				rr[0] = 1 / tms.rij(theta0, phi0p, theta1p, phi1);
				rr[1] = 1 / tms.rij(theta0, phi0n, theta1p, phi1);
				rr[2] = 1 / tms.rij(theta0, phi0p, theta1n, phi1);
				rr[3] = 1 / tms.rij(theta0, phi0n, theta1n, phi1);
				dg = r4h * (rr[0] - rr[1] - rr[2] + rr[3]);
				tst(2 * cy1, 2 * cy0 + 1) += dg;
				tst(2 * cy1, 2 * cy1 + 1) -= dg;
				tst(2 * cy0 + 1, 2 * cy1) += dg;
				tst(2 * cy1 + 1, 2 * cy1) -= dg;

				rr[0] = 1 / tms.rij(theta0, phi0p, theta1, phi1);
				rr[1] = 1 / tms.rij(theta0, phi0n, theta1, phi1);
				dg = rh * (rr[0] + rr[1] - 2 * r0);
				tst(2 * cy0 + 1, 2 * cy0 + 1) += dg;
				tst(2 * cy0 + 1, 2 * cy1 + 1) -= dg;
				tst(2 * cy1 + 1, 2 * cy0 + 1) -= dg;
				tst(2 * cy1 + 1, 2 * cy1 + 1) += dg;

			}
		}*/
	}
	void check()
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
		::printf("Shortest Distances: %lld, %.5e\n", mins.dim, mins.norm2() / mins.dim);
	}
	void run()
	{
		//tms.pos.printToTableTxt("E:\\files\\C++\\ComputePhysics\\A34\\Homework1_4\\pos.txt",true);
		for (unsigned long long c0(0); c0 < U.height; c0 += 2)
		{
			vec tp(U.data + (c0 + 1) * U.width4d, U.width4d, Type::Parasitic);
			double s(sin(tms.pos[c0]));
			tp /= s;
		}
		for (unsigned long long c0(0); c0 < U.width; c0 += 2)
		{
			double s(sin(tms.pos[c0]));
			for (unsigned long long c1(0); c1 < U.height; ++c1)
				U(c1, c0 + 1) /= s;
		}
		//U.print();
		//U.printToTableTxt("E:\\files\\C++\\ComputePhysics\\A34\\Homework1_4\\mat.txt");
		//must use L^TL
		mat tridiag(U.tridiagonalizationHouseholder());
		vec f(U.width, false);
		tridiag.implicitSymmetricQR(1e-42, f);
		fs = 1;
		freqs[0] = f[0];
		fn[0] = 1;
		for (unsigned long long c0(1); c0 < f.dim; ++c0)
		{
			unsigned long long c1(0);
			for (; c1 < fs; ++c1)
			{
				if (abs(f[c0] - freqs[c1] / fn[c1]) < 1e-8)
				{
					fn[c1]++;
					freqs[c1] += f[c0];
					break;
				}
			}
			if (c1 == fs)
			{
				freqs[fs] = f[c0];
				fn[fs++] = 1;
			}
		}
		freqs.dim = fs;
		fn.dim = fs;
		fs2r.dim = fs;
		freqs /= fn;
		for (unsigned long long c0(0); c0 < f.dim; ++c0)
			for (unsigned long long c1(0);; c1++)
			{
				double dt(f[c0] - freqs[c1]);
				if (abs(dt) < 1e-8)
				{
					fs2r[c1] += dt * dt;
					break;
				}
			}
		fs2r /= fn;
		fs2r.sqrt();
		freqs.abs().sqrt();
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
		tms.run(1e-20, rd, mt);
		::printf("%llu:\t%.11f\t%.3e\n", c0, tms.answer, tms.answer - tms.answers[c0]);
	}
	timer.end();
	timer.print();*/

	timer.begin();
	Vibration vbr(12, rd, mt);
	vbr.check();
	vbr.run();
	timer.end();
	timer.print();
	vbr.freqs.print(true);
	vbr.fn.print(true);
	vbr.fs2r.print(true);
}