#ifndef _M_MATRIX
#define _M_MATRIX
/*******************向量类*************************/


class Vector
{
public:
	Vector();				   //默认构造函数，向量大小为零
	Vector(int n);             //构造函数，指定向量大小
	Vector(const Vector &pVec);//拷贝构造函数
	~Vector();		           //析构函数
	
	int GetSize() const;     //获得向量大小

	double Norm(int i) const;//i模
	double Max_Norm() const; //无穷模

	const double &operator[](int n) const;     //对[]进行重载，以返回向量第n个元素,不能改变数据
	double &operator[](int n);                 //对[]进行重载，以返回向量第n个元素，可以改变数据
	double operator*(const Vector &pVec) const;//对*进行重载，以实现向量点乘
	
	Vector &operator=(const Vector &pVer); //对=进行重载，以实现向量赋值，大小为零的向量可以被任意大小的向量赋值
	
	void OutputData();		   //屏幕打印数据
	void OutputData(const char*Path,const char* op="w");//文件打印数据

	friend void V_Add(const Vector& left,const Vector& right,Vector& result);//向量加法
	friend void V_Sub(const Vector& left,const Vector& right,Vector& result);//向量减法
	friend void V_MulNum(double left,const Vector& right,Vector& result);//向量数乘 
	friend void V_MulNum(const Vector& left,double right,Vector& result);//向量数乘
protected:
	double* pVector;//向量指针
	int size;       //向量大小
};


/**********************矩阵类************************/

//矩阵元素结构体
typedef struct M_Element
{
	int pos;//元素列坐标(列序数)
	double value;//值
	struct M_Element *pNext;//指向下一个元素的指针
}M_Element;


class Matrix
{
public:
	
	Matrix();
	Matrix(int n);
	Matrix(const Matrix& matri);
	~Matrix();

	void SetVal(int m,int n,double value);//设置值
	void AddVal(int m,int n,double value);//追加值
	double GetVal(int m,int n) const;//取得值
	int GetSize()const{return size;};//获得矩阵大小
	int GetRealSize() const;//获得真实存储数据大小

	Matrix &operator=(const Matrix& matri);

	void OutputData();		   //屏幕打印数据
	void OutputData(char*Path);//文件打印数据
	void OutputRealData(char*Path);//文件打印真实存储数据
	
	void Grad(const Vector& right,Vector& X,double eps);//X传入初始值，传出计算结果
//	void Grad(const Vector& right,VEMFunction& X,double eps);
	void Pcg(const Vector& right,Vector& X,double eps);
	void Ldl(const Vector& right,Vector& X,double eps);
	void Ldl2(const Vector &right,Vector &X,double eps);
	void Ldl3(const Vector &right,Vector &X,double eps);

	friend void M_Add(const Matrix& left,const Matrix& right,Matrix& result);//矩阵加法,result必须事先指定正确的大小
	friend void M_Sub(const Matrix& left,const Matrix& right,Matrix& result);//矩阵减法,result必须事先指定正确的大小
	friend void M_MulNum(double left,const Matrix& right,Matrix& result);//矩阵数乘，矩阵为右操作数
	friend void M_MulNum(const Matrix& left,double right,Matrix& result);//矩阵数乘，矩阵为左操作数
	friend void M_MulVec(const Matrix& left,const Vector& right,Vector& result);//矩阵乘以向量

	int size;//矩阵大小，size*size方阵
	M_Element* Header;//每行头指针数组
};
#endif
