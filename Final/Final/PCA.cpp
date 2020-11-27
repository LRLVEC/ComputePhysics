#include <_BLAS.h>
#include <_Time.h>
#include <_File.h>
#include <random>

using namespace BLAS;

void readData(char const* str, mat& a, mat& aT)
{
	int n(0), dn;
	for (unsigned long long c0(0); c0 < 10000; ++c0)
		for (unsigned long long c1(0); c1 < 6; ++c1)
		{
			sscanf(str + n, "%lf%n", &aT(c1, c0), &dn);
			n += dn;
		}
	for (unsigned long long c0(0); c0 < 6; ++c0)
	{
		vec v(aT.row(c0));
		v -= v.average();
	}
	for (unsigned long long c0(0); c0 < 10000; ++c0)
		for (unsigned long long c1(0); c1 < 6; ++c1)
			a(c0, c1) = aT(c1, c0);
}
void writeAnswer(mat& const answer, File& file)
{
	String<char>tmp;
	tmp.malloc(600000);
	char t[128];
	for (unsigned long long c0(0); c0 < 10000; ++c0)
	{
		sprintf(t, "%.4f %.4f %.4f %.4f %.4f %.4f\n",
			answer(0, c0), answer(1, c0), answer(2, c0),
			answer(3, c0), answer(4, c0), answer(5, c0));
		tmp += t;
	}
	file.createText("answer.txt", tmp);
}

int main()
{
	std::mt19937 mt(time(nullptr));
	std::uniform_real_distribution<double> rd(0, 1.0);
	std::uniform_int_distribution<unsigned long long> rduint(1, 10);
	Timer timer;
	File file("./");
	String<char>str(file.readText("data.txt"));

	mat A(6, 10000, false);
	mat AT(10000, 6, false);
	mat ATA(6, 6, false);
	//mat AAT(10000, 10000, false);
	vec eigenvalues(6, false);
	mat eigenvectorsATA(6, 6, false);
	mat eigenvectorsAAT(10000, 6, false);
	mat answers(6, 10000, false);
	vec variance(6, false);
	timer.begin();
	readData(str, A, AT);
	timer.end();
	timer.print("Read data: ");

	timer.begin();
	AT(A, ATA);
	timer.end();
	timer.print("AT*A: ");

	//timer.begin();
	//A(AT, AAT);
	//timer.end();
	//timer.print("A*AT: ");

	ATA /= 10000;
	mat ATABackup(ATA);
	ATA.print();
	ATA.printToTableTxt("./a.txt");

	mat triATA(ATABackup.tridiagonalizationHouseholder());
	triATA.implicitSymmetricQR(1e-70, eigenvalues);
	eigenvalues.qsortD().print();
	eigenvalues.printToTableTxt("eigenvalues.txt", false);

	timer.begin();
	ATA.inversePowerEigenvectors(eigenvalues, eigenvectorsATA);
	timer.end();
	timer.print("EigenvectorsATA: ");
	eigenvectorsATA.printToTableTxt("./eigenvectorsATA.txt");
	eigenvectorsATA.print();

	mat answer(eigenvectorsATA(AT));
	writeAnswer(answer, file);

	

	//vec tst(6, false);
	//mat tmt(6, 6, false);
	//ATA.powerMaxEigenvector(tst);
	//tst.print();
	//tmt.row(0) = tst;
	//ATA.powerEigenvectors(tmt);
	//tmt.print();

	//timer.begin();
	//AAT.inversePower(eigenValues, eigenVectorsAAT);
	//timer.end();
	//timer.print("EigenVectorsAAT: ");
	//eigenVectorsAAT.print();

	//mat test(2048, 2048, false);
	//randomMatSymmetric(test, mt, rd, 0.0002);
	//timer.begin();
	//test.solveCholesky();
	//timer.end();
	//timer.print();
}