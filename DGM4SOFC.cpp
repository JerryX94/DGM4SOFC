// DGM4SOFC_v4.2.cpp: 主项目文件。
// 在v4.1的基础上，引入初始可逆电势计算，增加了备注。

#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#define INPUTFILENAME1 "inputie.dat"
#define INPUTFILENAME2 "input.dat"
#define OUTPUTFILENAME "output.dat"
#define ALPHA 1.//Anodic transfer coefficient
#define BETA .5//Cathodic transfer coefficient
#define D 2.E-6//m
#define EACT 1.E5//J/mol
#define ERRMAX 1.E-6
#define F 96485.//C/mol
#define K1 1.92//Anodic Thickness Coefficient
#define K2 1.3//Cathodic Thickness Coefficient
#define KP .01;//WGSR Coefficient
#define MAXIT 10000//Max Iterature Num.
#define MAXNIE 30//Operation Current Density Input Num.
#define N 2//Component Num.
#define NDISC 100//Discretization Num.
#define NREACT 1//Eletrochemical Reaction Num.
#define O20 .7//Oxygen Molar Fraction in Cathode Channel
#define PA 1.013E5//Pa
#define PA0 1.013E5//Pa
#define PI 3.1416
#define PORO .45
#define R 8.314//J/mol.K
#define T 1073.15//K
#define TREF 1073.15//K
#define THICKA 1.1E-3//m
#define THICKE 1.E-5//m
#define PODIVTO .1

int InIe(double *IeQueue);
void Init();
void InWg();
double CBg();
double CDKn(int species);
double CDij(int speciesi, int speciesj);
int Iter(double Bg, double DKn[], double Dij[][N]);
void Exha(int speciesi);
void CWg(double Z);
double CMu(double *Xave);
double Err(double *X0, double *X1, double *DX, double P0, double P1, int it);
double Ei0(int ri);
double CEta();
double EtaCathode();
void Enqueue(int it, int iIE, double OutputQueue[][N + 4]);
void Outp(int nIE, double OutputQueue[][N + 4]);

struct Gas{
public:
	double P0;//Pa
	double DP;//Pa
	double DPlist[NDISC];//Pa
	double X0[N];
	double X1[N];
	double Xlist[N][NDISC];
	double TotalNg;//mol/m2.s
	double Ng[N];//mol/m2.s
	double DNg[N];//mol/m2.s
	double Mg[N];//g/mol
	double Vd[N];//cm3/mol
	double ReactOrder[N];
	double J0ref[N];//A/m3
	double ViscCoef[N][7];
};

struct React1{
public:
	short Ne;
	short Reactant;
	short Product;
	double w;//Species Current Factor
	double g[3];
	double s[3];
};

struct React2{
public:
	short Switch;
	short Reactant[2];
	short Product[2];
	double Keq;
};

Gas mixture;
React1 elechem[NREACT];
React2 wgsr;
double IE = 0.;//电流密度，A/m2；
double MAXOVERP = 0.;//SOFC可逆电压，V；

int main(){
	//程序主函数；
	double IeQueue[MAXNIE];
	double OutputQueue[MAXNIE][N + 4];
	double DKn[N];//Knudsen扩散系数，m2/s
	double Dij[N][N] = { 0 };//二元（连续）扩散系数，m2/s
	int nIE = InIe(IeQueue);
	for (int iIE = 0; iIE < nIE; iIE++){
		IE = IeQueue[iIE];
		Init();
		double Bg = CBg();
		for (int i = 0; i < N; i++)
			DKn[i] = CDKn(i);
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++){
				if (i != j) Dij[i][j] = Dij[j][i] = CDij(i, j);
			}
		Enqueue(Iter(Bg, DKn, Dij), iIE, OutputQueue);
	}
	Outp(nIE, OutputQueue);
	return 0;
}

int InIe(double *IeQueue){
	//本函数从文件读取一组电流密度输入，输入数据单位为A/cm2；
	using namespace std;
	ifstream fin(INPUTFILENAME1);
	int nIE;
	fin >> nIE;
	for (int i = 0; i < nIE; i++){
		fin >> IeQueue[i];
		IeQueue[i] *= 10000.;//单位变换，A/cm2->A/m2；
	}
	fin.close();
	return nIE;
}

void Init(){
	//本函数读取文件流输入，对计算所需的各变量进行初始化；
	using namespace std;
	ifstream fin(INPUTFILENAME2);
	for (int i = 0; i < NREACT; i++){
		fin >> elechem[i].Ne >> elechem[i].Reactant >> elechem[i].Product;
	}
	fin >> wgsr.Switch;
	if (wgsr.Switch){
		fin >> wgsr.Reactant[0] >> wgsr.Reactant[1];
		fin >> wgsr.Product[0] >> wgsr.Product[1];
		InWg();
	}
	mixture.P0 = PA;
	mixture.DP = 0.;
	for (int i = 0; i < N; i++){
		fin >> mixture.X0[i];
		mixture.X1[i] = mixture.X0[i];
	}
	for (int i = 0; i < N; i++){
		mixture.DNg[i] = 0.;
		fin >> mixture.Ng[i];//Input Stoichiometric Coefficient
		short flag = 0;
		for (int j = 0; j < NREACT; j++){
			if ((elechem[j].Reactant == i) || (elechem[j].Product == i)){
				mixture.Ng[i] *= IE / (elechem[j].Ne * F);
				flag++;
				break;
			}
		}
		if (!flag) mixture.Ng[i] = 0.;
	}
	for (int i = 0; i < N; fin >> mixture.Mg[i++]);
	for (int i = 0; i < N; fin >> mixture.Vd[i++]);
	for (int i = 0; i < N; fin >> mixture.ReactOrder[i++]);
	for (int i = 0; i < N; fin >> mixture.J0ref[i++]);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < 7; j++)
			fin >> mixture.ViscCoef[i][j];
	double rt = 0.;
	for (int i = 0; i < NREACT; i++){
		rt += mixture.X0[elechem[i].Reactant] * elechem[i].Ne;
	}
	mixture.TotalNg = 0.;
	for (int i = 0; i < NREACT; i++){
		elechem[i].w = mixture.X0[elechem[i].Reactant] * elechem[i].Ne / rt;
		mixture.Ng[elechem[i].Reactant] *= elechem[i].w;
		mixture.Ng[elechem[i].Product] *= elechem[i].w;
		mixture.TotalNg += mixture.Ng[elechem[i].Reactant] * elechem[i].Ne;
	}
	for (int i = 0; i < NREACT; i++){
		for (int j = 0; j < 3; j++){
			fin >> elechem[i].g[j];
		}
		for (int j = 0; j < 3; j++){
			fin >> elechem[i].s[j];
		}
	}
	fin.close();
	return;
}

void InWg(){
	//本函数对水气转换反应（WGSR）进行初始化，计算反应化学平衡常数；
	double Z = 1000. / T - 1.;
	wgsr.Keq = -.2935 * Z * Z * Z + .635 * Z * Z + 4.1788 * Z + .3169;
	wgsr.Keq = exp(wgsr.Keq);
}

double CBg(){
	//本函数计算多孔介质渗透率；
	//不同文献中渗透率计算公式有一定差异，但对计算结果影响有限；
	double Bg = (pow(D, 2) / 180.) / (PORO * PORO / pow(1 - PORO, 2));
	//double Bg = 1.7E-10;
	return Bg;
}

double CDKn(int species){
	//本函数计算一种气体组分的Knudsen扩散系数；
	double eDKn = D / 3. * sqrt(8. * R * T\
		/ (PI * mixture.Mg[species] / 1000.));
	eDKn *= PODIVTO;
	return eDKn;
}

double CDij(int speciesi, int speciesj){
	//本函数计算气体i与气体j之间的二元（连续）扩散系数；
	double k = 1.0E-3;
	double mij = 1. / mixture.Mg[speciesi] + 1. / mixture.Mg[speciesj];
	double nume = k * pow(T, 1.75) * sqrt(mij);
	double deno = pow(mixture.Vd[speciesi], 1. / 3)\
		+ pow(mixture.Vd[speciesj], 1. / 3);
	deno = deno * deno * PA / PA0;
	double eDij = PODIVTO * nume / deno * 1.E-4;
	return eDij;
}

int Iter(double Bg, double DKn[], double Dij[][N]){
	//本函数基于DGM（尘气模型）计算SOFC阳极质量传输过程；
	int it = 0;//迭代次数；
	const double Z = K1 * THICKA / NDISC;//离散等效电极厚度，凑实验数据；
	for (int k = 0; k < NDISC; k++){
		double DX[N] = { 0 };//摩尔因数差；
		double Xave[N];//平均摩尔因数；
		double Ngave[N];//平均流量；
		double DPnew = mixture.DP / NDISC;//更新压力差；
		double DPold = DPnew;
		double Pave;//平均压力；
		do{
			if (it > MAXIT){
				it = -1;
				break;
			}
			//水气转换反应，先迭代两次使浓度分布偏差减小后再行实施；
			if (wgsr.Switch){
				if (it > 1) CWg(Z);
			}
			//调整权重系数；
			for (int i = 0; i < NREACT; i++){
				elechem[i].w = mixture.Ng[elechem[i].Reactant] * elechem[i].Ne\
					/ mixture.TotalNg;
			}
			//获取新值；
			for (int i = 0; i < N; i++){
				mixture.X1[i] = mixture.X0[i] + DX[i];
				if (mixture.X1[i] < 0.){
					mixture.X1[0] += mixture.X1[i];
					mixture.X1[i] = 0.;
					Exha(i);
				}
			}
			//调整权重系数；
			for (int i = 0; i < NREACT; i++){
				elechem[i].w = mixture.Ng[elechem[i].Reactant] * elechem[i].Ne\
					/ mixture.TotalNg;
			}
			//更新计算指标；
			DPold = DPnew;
			for (int i = 0; i < N; i++){
				Xave[i] = .5 * (mixture.X0[i] + mixture.X1[i]);
				Ngave[i] = mixture.Ng[i] - .5 * mixture.DNg[i];
			}
			Pave = mixture.P0 + .5 * DPold;
			double Mu = CMu(Xave);//重新计算粘度；
			//计算压力差；
			double nume = 0.;
			double deno = 1.;
			for (int i = 0; i < N; i++)
				nume += Ngave[i] / DKn[i];
			for (int i = 0; i < N; i++)
				deno += (Bg / Mu) * Pave * Xave[i] / DKn[i];
			DPnew = -(nume / deno) * Z * R * T;
			//计算浓度差；
			for (int i = 1; i < N; i++){
				DX[i] = 0.;
				for (int j = 0; j < N; j++){
					if (j != i){
						nume = Xave[j] * Ngave[i] - Xave[i] * Ngave[j];
						deno = Dij[i][j];
						DX[i] += nume / deno;
					}
				}
				DX[i] = R * T * (DX[i] + Ngave[i] / DKn[i]) + Xave[i]\
					* DPold * (1. + Bg * Pave / (Mu * DKn[i])) / Z;
				DX[i] *= -Z / Pave;
			}
			DX[0] = 0.;
			for (int i = 1; i < N; DX[0] -= DX[i++]);
		} while (Err(mixture.X0, mixture.X1, DX, DPold, DPnew, ++it)\
			> ERRMAX);
		for (int i = 0; i < N; i++){
			mixture.Xlist[i][k] = mixture.X0[i];
			mixture.X0[i] = mixture.X1[i];
		}
		mixture.DPlist[k] = DPnew;
		mixture.DP += DPnew;
	}
	for (int i = 0; i < N; i++){
		mixture.X0[i] = mixture.Xlist[i][0];
	}
	return it;
}

double Err(double *X0, double *X1, double *DX, double DP0, double DP1, int it){
	//本函数计算DGM迭代过程中各组分摩尔分数相对上一轮迭代值的偏差；
	double Err2 = 0.;
	for (int i = 0; i < N; i++){
		double diff = X0[i] + DX[i] - X1[i];
		Err2 += (diff * diff) / (X0[i] * X0[i] + 1.E-10);
	}
	double diff = DP0 - DP1;
	Err2 += (diff * diff) / (DP0 * DP0 + 1.E-10);
	if (!(it % 50))
		std::cout << "Iteration = " << it << ", Err = " << Err2 << "\n";
	return Err2;
}

void Exha(int speciesi){
	double relax = .99;//松弛因子，修正流量；
	double dNg = 0.;
	for (int i = 0; i < NREACT; i++)
		if (elechem[i].Reactant == speciesi){
			dNg += mixture.Ng[elechem[i].Reactant] * elechem[i].Ne\
				/ elechem[0].Ne * (1. - relax);
			mixture.Ng[elechem[i].Reactant] *= relax;
			mixture.Ng[elechem[i].Product] *= relax;
		}
	//用化学反应elechem[0]补偿因修正而减少的总流量；
	//因此必须始终保证化学反应0的反应物最后一个耗尽；
	mixture.Ng[elechem[0].Reactant] += dNg;
	mixture.Ng[elechem[0].Product] -= dNg;
	return;
}

void CWg(double Z){//略有问题，容易导致发散
	const int WGSplit = 5;//空间离散节点数
	double r[WGSplit];
	for (int i = 0; i < WGSplit; i++){
		double x1 = mixture.X0[wgsr.Product[0]]\
			+ i * (mixture.X1[wgsr.Product[0]]\
			- mixture.X0[wgsr.Product[0]]) / (WGSplit - 1);
		double x2 = mixture.X0[wgsr.Reactant[0]]\
			+ i * (mixture.X1[wgsr.Reactant[0]]\
			- mixture.X0[wgsr.Reactant[0]]) / (WGSplit - 1);
		double x3 = mixture.X0[wgsr.Reactant[1]]\
			+ i * (mixture.X1[wgsr.Reactant[1]]\
			- mixture.X0[wgsr.Reactant[1]]) / (WGSplit - 1);
		double x4 = mixture.X0[wgsr.Product[1]]\
			+ i * (mixture.X1[wgsr.Product[1]]\
			- mixture.X0[wgsr.Product[1]]) / (WGSplit - 1);
		r[i] = .0171 * exp(-103191. / (R * T))\
			* (x2 * x3 - x1 * x4 / wgsr.Keq) * PA * PA;
	}
	double rtotal = 0.;
	for (int i = 0; i < WGSplit; rtotal += r[i++]);
	rtotal *= KP;
	for (int i = 0; i < 2; i++){
		mixture.DNg[wgsr.Reactant[i]] = -rtotal / WGSplit * Z;
		mixture.DNg[wgsr.Product[i]] = rtotal / WGSplit * Z;
	}
	return;
}

double CMu(double *Xave){
	//本函数计算气体混合物动力粘性系数；
	double MuSp[N] = { 0 };
	double phi[N][N] = { 0 };
	for (int i = 0; i < N; i++)
		for (int j = 0; j < 7; j++)
			MuSp[i] += mixture.ViscCoef[i][j] * pow(T / 1000., j);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			if (j != i){
				double nume = 1.;
				nume += sqrt(MuSp[i] / MuSp[j]) * pow(mixture.Mg[j]\
					/ mixture.Mg[i], .25);
				nume *= nume;
				double deno = sqrt(8. * (1 + mixture.Mg[i]\
					/ mixture.Mg[j]));
				phi[i][j] = nume / deno;
			}
	double Mum = 0.;
	for (int i = 0; i < N; i++){
		double deno = 1.;
		for (int j = 0; j < N; j++){
			if (j != i) deno += Xave[j] * phi[i][j] / Xave[i];
		}
		Mum += MuSp[i] / deno * 1.E-7;//单位变换：微泊(muP)->Pa.s；
	}
	return Mum;
}

double Ei0(int ri){
	//本函数计算电化学反应elechem[ri]提供的可逆电压；
	const double Eff = .95;//燃料电池效率，根据实验值对比修正；
	const double T0 = 1000.;
	double E0 = -(elechem[ri].g[1] - elechem[ri].g[0]\
		- .5 * elechem[ri].g[2]) * 1.E3 / (elechem[ri].Ne * F);
	double Et = E0 + (elechem[ri].s[1] - elechem[ri].s[0]\
		- .5 * elechem[ri].s[2]) / (elechem[ri].Ne * F) * (T - T0);
	double E = Eff * (Et - (R * T) / (elechem[ri].Ne * F)\
		* log(mixture.X0[elechem[ri].Product]\
		/ mixture.X0[elechem[ri].Reactant] / sqrt(O20)));
	return E;
}

double CEta(){
	//本函数通过DGM解得的TPB中各组分浓度变化，计算SOFC极化损失及欧姆损失；
	//由于本函数依赖的表达式[Bao, 2016]是隐式的，因此求解方法为迭代；
	double TotalEta = 1.E10;
	double SpeciesEta = 0.;
	double OhmEta = 0.;
	double IonCond = 3.34E4 * exp(-10300. / T);
	MAXOVERP = 0.;
	for (int i = 0; i < NREACT; i++){
		MAXOVERP += Ei0(i) * elechem[i].w;
	}
	for (int i = 0; i < NREACT; i++){
		double I0refb = IonCond * R * T\
			/ ((ALPHA + BETA) * elechem[i].Ne * F);
		I0refb *= mixture.J0ref[elechem[i].Reactant];
		I0refb *= pow(PA * mixture.X0[elechem[i].Reactant] / PA0, \
			mixture.ReactOrder[elechem[i].Reactant]);
		I0refb *= pow(PA * mixture.X0[elechem[i].Product] / PA0, \
			mixture.ReactOrder[elechem[i].Product]);
		I0refb = sqrt(I0refb);
		double I0 = I0refb * exp(-EACT / R * (1. / T - 1 / TREF));
		double Icomp = elechem[i].w * IE;//各电化学反应贡献的电流；
		double sub0, sub = 0.;
		do{
			sub0 = sub;
			sub = exp(-BETA * sub) * mixture.X1[elechem[i].Product]\
				/ mixture.X0[elechem[i].Product] + Icomp\
				/ (I0 * mixture.X0[elechem[i].Reactant]);
			sub = log(sub * mixture.X0[elechem[i].Reactant]\
				/ mixture.X1[elechem[i].Reactant]) / ALPHA;
		} while (abs(sub - sub0) / (sub0 + 1.E-10) > ERRMAX);
		SpeciesEta = sub * R * T / (elechem[i].Ne * F);
		if (SpeciesEta > MAXOVERP) SpeciesEta = MAXOVERP;
		if (TotalEta > SpeciesEta) TotalEta = SpeciesEta;//取最小，可能有问题；
	}
	OhmEta = IE * THICKE / IonCond;//欧姆损失；
	TotalEta += OhmEta;
	return TotalEta;
}

double EtaCathode(){
	//本函数与阳极计算类似，通过DGM计算阴极质量传输和极化损失
	struct CathodeGas{
		double X0;
		double X1;
		double Ng;//mol/m2.s
		double Mg;//g/mol
		double Vd;//cm3/mol
		double EDKn;
		double MuSp;
		double ViscCoef[7];
	};
	struct CathodeGas O2 = {
		O20, O20, .5 * IE / (2. * F), 31.999, 16.3, 0., 0.,
		{ -1.6918, 889.75, -892.79, 905.98, -598.36, 221.64, -34.754 }
	};
	struct CathodeGas N2 = {
		1. - O20, 1. - O20, 0., 28.013, 18.5, 0., 0.,
		{ 1.2719, 771.45, -809.20, 832.47, -553.93, 206.15, -32.430 }
	};

	O2.EDKn = PODIVTO * (D / 3. * sqrt(8. * R * T / (PI * O2.Mg / 1000.)));
	N2.EDKn = PODIVTO * (D / 3. * sqrt(8. * R * T / (PI * N2.Mg / 1000.)));
	const double Z = K2 * THICKA;
	double nume = 1. / O2.Mg + 1. / N2.Mg;
	nume = 1.E-3 * pow(T, 1.75) * sqrt(nume);
	double deno = pow(O2.Vd, 1. / 3) + pow(N2.Vd, 1. / 3);
	deno = deno * deno * PA / PA0;
	double EDij = PODIVTO * nume / deno * 1.E-4;
	//DGM计算阴极质量传输
	double DX[2] = { 0 };
	double Xave[2];
	double Bg = CBg();
	double P0 = PA;
	double DPold = 0.;
	double DPnew = 0.;
	double Pave, Err2;
	do{
		//获取新值
		O2.X1 = O2.X0 + DX[0];
		N2.X1 = N2.X0 + DX[1];
		DPold = DPnew;
		Xave[0] = .5 * (O2.X0 + O2.X1);
		Xave[1] = .5 * (N2.X0 + N2.X1);
		Pave = P0 + .5 * DPold;
		//重新计算粘度
		O2.MuSp = N2.MuSp = 0.;
		for (int i = 0; i < 7; i++){
			O2.MuSp += O2.ViscCoef[i] * pow(T / 1000., i);
			N2.MuSp += N2.ViscCoef[i] * pow(T / 1000., i);
		}
		nume = 1. + sqrt(O2.MuSp / N2.MuSp) * pow(N2.Mg / O2.Mg, .25);
		nume *= nume;
		deno = sqrt(8. * (1 + O2.Mg / N2.Mg));
		double phiij = nume / deno;
		nume = 1. + sqrt(N2.MuSp / O2.MuSp) * pow(O2.Mg / N2.Mg, .25);
		nume *= nume;
		deno = sqrt(8. * (1 + N2.Mg / O2.Mg));
		double phiji = nume / deno;
		double Mu = O2.MuSp / (1. + Xave[1] * phiij / Xave[0])\
			+ N2.MuSp / (1. + Xave[0] * phiji / Xave[1]);
		Mu *= 1.E-7;//单位变换：微泊(muP)->Pa.s
		//计算压力差
		nume = O2.Ng / O2.EDKn + O2.Ng / O2.EDKn;
		deno = 1. + (Bg / Mu) * Pave * (1. / O2.EDKn + 1. / N2.EDKn);
		DPnew = -(nume / deno) * Z * R * T;
		//计算浓度差
		DX[0] = (Xave[1] * O2.Ng - Xave[0] * N2.Ng) / EDij;
		DX[0] = -Z / Pave * (R * T * (DX[0] + O2.Ng / O2.EDKn) + Xave[0]\
			* DPold * (1. + Bg * Pave / (Mu * O2.EDKn)) / Z);
		DX[1] = -DX[0];
		//计算误差
		Err2 = 0.;
		double diff = O2.X0 + DX[0] - O2.X1;
		Err2 += (diff * diff) / (O2.X0 * O2.X0 + 1.E-10);
		diff = N2.X0 + DX[1] - N2.X1;
		Err2 += (diff * diff) / (N2.X0 * N2.X0 + 1.E-10);
		diff = DPold - DPnew;
		Err2 += (diff * diff) / (DPold * DPold + 1.E-10);
	} while (Err2 > ERRMAX);
	//计算阴极极化损失
	double sub0, sub = 0.;
	do{
		sub0 = sub;
		sub = exp(-sub) + IE / (6945.7 * pow(Pave * O2.X0 / PA, .25));
		sub = log(sub * O2.X0 / O2.X1);
	} while (abs(sub - sub0) / (sub0 + 1.E-10) > ERRMAX);
	double Eta = sub * R * T / F;
	if (Eta > MAXOVERP) Eta = MAXOVERP;
	return Eta;
}

void Enqueue(int it, int iIE, double OutputQueue[][N + 4]){
	//本函数通过队列存储多个电流密度下的多组计算数据；
	double Eta = CEta() + EtaCathode();
	OutputQueue[iIE][0] = (double)it;
	for (int i = 0; i < N; i++)
		OutputQueue[iIE][i + 1] = mixture.X1[i];
	OutputQueue[iIE][N + 1] = mixture.DP;
	OutputQueue[iIE][N + 2] = Eta;
	OutputQueue[iIE][N + 3] = MAXOVERP - Eta;
}

void Outp(int nIE, double OutputQueue[][N + 4]){
	//本函数对程序计算结果进行格式化文件流输出；
	using namespace std;
	ofstream fout(OUTPUTFILENAME);
	fout << "Iter_Times ";
	for (int i = 0; i < N; i++){
		fout << "AE_Interf_X" << i << " ";
	}
	fout << "Pressure_Dif " << "Concen_Overp " << "Voltage" << endl;
	fout.setf(ios::left);
	for (int k = 0; k < nIE; k++){
		fout << setw(11) << (int)OutputQueue[k][0];
		for (int i = 1; i < N + 4; fout << setw(13) << OutputQueue[k][i++]);
		fout << endl;
	}
	fout.close();
	return;
}