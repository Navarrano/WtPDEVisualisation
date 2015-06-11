#ifndef PTR
#define PTR


template<class T> T* clone(const T* tp)
{
	return tp->clone();
}


template<class T> 
class Ptr {

private:
	T* p;
	size_t* refptr;
	static int count;
public:
	Ptr();
	Ptr(T* t);
	Ptr(const Ptr<T>& h);
	Ptr<T>& operator=(const Ptr<T>& rhs);
	~Ptr();
	operator bool() const;
	T& operator*() const;
	T* operator->() const;
	void make_unique();
	static int getCount() { return count; }
};


template<class T>
Ptr<T>::Ptr() : refptr(new size_t(1)), p(0) { count++; } /* std::cout << "default object created\n"; std::cout << count << "\n";*/ //cn++; std::cout << cn << "\n";}

template<class T>
Ptr<T>::Ptr(T* t) : refptr(new size_t(1)), p(t) { count++; }/* std::cout << "object created\n"; std::cout << count << "\n" */; //cn++; std::cout << cn << "\n";}

template<class T>
Ptr<T>::Ptr(const Ptr<T>& h) : refptr(h.refptr), p(h.p) { ++*refptr; count++; }/* std::cout << "object copied copy constr\n"; std::cout << count << "\n"; */ //cn++; std::cout << cn << "\n";}

template<class T>
Ptr<T>& Ptr<T>::operator=(const Ptr<T>& rhs) 
{
	++*rhs.refptr;
	// free left-hand side, destroying pointers if appropriate
	if (--*refptr == 0) {
		delete refptr;
		delete p;
	}
	// copy in values from the right-hand side
	refptr = rhs.refptr;
	p = rhs.p;
	return *this;
}

template<class T>
Ptr<T>::~Ptr()
{
	if (--*refptr == 0) {
		delete refptr;
		delete p;
	}
	count--;
//	cn--;
//	std::cerr << "object destroyed\n";
//	std::cout << cn << "\n";
}


template<class T>
Ptr<T>::operator bool() const
{
	return p!=0;
}

template<class T>
T& Ptr<T>::operator*() const
{
	if (p)
		return *p;
	throw std::runtime_error("unbound Ref_handle");
}

template<class T>
T* Ptr<T>::operator->() const
{
	if (p)
		return p;
	throw std::runtime_error("unbound Ref_handle");
}

template<class T>
void Ptr<T>::make_unique()
{
	if (*refptr != 1) {
		--refptr;
		refptr = new size_t(1);
		p = p ? clone() : 0;
	}
}

template<class T>
int Ptr<T>::count =  0;


#endif