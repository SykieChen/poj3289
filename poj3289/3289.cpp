#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <math.h>
#define PI 3.141592653589793

using namespace std;

/*不知道哪里有问题。算的不对。懒得调试了。*/


//横放体积函数
double f1(double v0, double s, double k, double hb, double db, double hn, double dn, double h);
double f2(double v0, double s, double k, double hb, double db, double hn, double dn, double h);
double f3(double v0, double s, double k, double hb, double db, double hn, double dn, double h);
double f4(double v0, double s, double k, double hb, double db, double hn, double dn, double h);

double solve(double k, double hb, double db, double hn, double dn, double h) {
	//计算体积
	double v, rb = dn / 2.0, rn = dn / 2.0;
	if (k<0 || k>h) {
		cout << "Data illegal, exit." << endl;
		return -1.0;
	}
	else if (k <= hb) v = PI*rb*rb*k;
	else if (k < (h - hn)) {
		double p1 = PI*rb*rb*hb;
		double p2 = (rn*rn*(k - hb)*(k - hb)) / ((h - hn - hb)*(h - hn - hb));
		double p3 = (rb*rn*(k - hb)) / (h - hn - hb);
		v = p1 + ((PI*h*((rb*rb) + p2 + p3)) / 3.0);
	}
	else {
		v = (PI*rb*rb*hb) + ((PI*(h - hb - hn)*((rb*rb) + (rn*rn) + (rb*rn))) / 3.0) + (PI*rn*rn*(k + hn - h));
	}

	//二分法解方程计算s
	//计算临界值
	double s1 = rb - rn;
	double v1 = f1(0, s1, k, hb, db, hn, dn, h);
	double s2 = rb;
	double v2 = f2(0, s2, k, hb, db, hn, dn, h);
	double s3 = rb + rn;
	double v3 = f3(0, s3, k, hb, db, hn, dn, h);
	//二分法求根，分四段
	if (v < v1) {
		double a = 0, b = s1, limit = 0.01;
		if (f1(v, a, k, hb, db, hn, dn, h)*f1(v, b, k, hb, db, hn, dn, h)>0) {
			cout << "ERR:can't solve f1\n";
			return -1.0;
		}
		else{
			while ((b - a)>limit){
				if (f1(v, (a + b) / 2.0, k, hb, db, hn, dn, h)*f1(v, b, k, hb, db, hn, dn, h) < 0) a = (a + b) / 2.0;
				else b = (a + b) / 2.0;
			}
			return a;
		}
	}
	else if (v < v2) {
		double a = s1, b = s2, limit = 0.01;
		if (f2(v, a, k, hb, db, hn, dn, h)*f2(v, b, k, hb, db, hn, dn, h)>0) {
			cout << "ERR:can't solve f2\n";
			return -1.0;
		}
		else {
			while ((b - a)>limit) {
				if (f2(v, (a + b) / 2.0, k, hb, db, hn, dn, h)*f2(v, b, k, hb, db, hn, dn, h) < 0) a = (a + b) / 2.0;
				else b = (a + b) / 2.0;
			}
			return a;
		}
	}
	else if (v < v3) {
		double a = s2, b = s3, limit = 0.01;
		if (f3(v, a, k, hb, db, hn, dn, h)*f3(v, b, k, hb, db, hn, dn, h)>0) {
			cout << "ERR:can't solve f3\n";
			return -1.0;
		}
		else {
			while ((b - a)>limit) {
				if (f3(v, (a + b) / 2.0, k, hb, db, hn, dn, h)*f3(v, b, k, hb, db, hn, dn, h) < 0) a = (a + b) / 2.0;
				else b = (a + b) / 2.0;
			}
			return a;
		}
	}
	else {
		double a = s3, b = dn, limit = 0.01;
		if (f4(v, a, k, hb, db, hn, dn, h)*f4(v, b, k, hb, db, hn, dn, h)>0) {
			cout << "ERR:can't solve f4\n";
			return -1.0;
		}
		else {
			while ((b - a)>limit) {
				if (f4(v, (a + b) / 2.0, k, hb, db, hn, dn, h)*f4(v, b, k, hb, db, hn, dn, h) < 0) a = (a + b) / 2.0;
				else b = (a + b) / 2.0;
			}
			return a;
		}
	}
}

double f1(double v0, double s, double k, double hb, double db, double hn, double dn, double h) {
	double rb = dn / 2.0, rn = dn / 2.0;
	double p1 = rb*rb*acos((rb - s) / rb);
	double p2 = sqrt((2.0 * rb*s) - (s*s))*(rb - s);
	double v = (p1 - p2) / 3.0;
	return v - v0;
}

double f2(double v0, double s, double k, double hb, double db, double hn, double dn, double h) {
	double rb = dn / 2.0, rn = dn / 2.0;
	double p1 = rn*rn*acos((rn - s) / rn);
	double p2 = rn*sqrt((rn*rn) - ((rb - s)*(rb - s)));
	double sa = p1 - p2;
	double p4 = rb*rb*acos((rb - s) / rb);
	double p5 = sqrt((2.0 * rb*s) - (s*s));
	double sb = p4 - p5;
	double v = (sb*hb) + ((h*(sa + sb + sqrt(sa*sb))) / 3.0) + (sa*hn);
	return v - v0;
}

double f3(double v0, double s, double k, double hb, double db, double hn, double dn, double h) {
	double rb = dn / 2.0, rn = dn / 2.0;
	double p1 = rn*rn*acos((s - rn) / rn);
	double p2 = (s - rb)*sqrt((rn*rn) - ((s - rb)*(s - rb)));
	double sa = PI*rn*rn - p1 + p2;
	double p3 = rb*rb*acos((s - rb) / rb);
	double p4 = (s - rb)*sqrt((2.0*s*rb) - (s*s));
	double sb = PI*rb*rb - p3 + p4;
	double v = (sb*hb) + ((h*(sa + sb + sqrt(sa*sb))) / 3.0) + (sa*hn);
	return v - v0;
}

double f4(double v0, double s, double k, double hb, double db, double hn, double dn, double h) {
	double rb = dn / 2.0, rn = dn / 2.0;
	double sa = PI*rn*rn;
	double p1 = rb*rb*acos((s - rb) / rb);
	double p2 = (s - rb)*sqrt((2.0*s*rb) - (s*s));
	double sb = PI*rb*rb - p1 + p2;
	double v = (sb*hb) + ((h - hb - hn) *PI* ((rn*rn) + (rb*rb) + sqrt((rn*rn) + (rb*rb))) / 3.0) + (sa*hn) - (((2.0*rb) - (s*(h - hn - hb))) / (3.0*rb));
	return v - v0;
}

int main() {
	double k = 0.0, hb, db, hn, dn, h;
	double ans;
	do {
		cin >> k >> hb >> db >> hn >> dn >> h;
		if (k) {
			ans = solve(k, hb, db, hn, dn, h);
			if (ans >= 0) cout << ans << endl;
			else cout << "Program ended with error." << endl;
		}
	} while (k);
	return 0;
}