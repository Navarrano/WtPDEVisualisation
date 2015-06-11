#ifndef ABSTRACT_PANEL_H
#define ABSTRACT_PANEL_H

#include <Wt/WPanel>
#include <Wt/WContainerWidget>
#include <string>
#include <Wt/Chart/WCartesian3DChart>
#include <boost/signal.hpp>

#include "headers.h"

using namespace Wt;



typedef boost::signal <void(ChartData)> CalcSignal;

class AbstractPanel : public Wt::WPanel
{
public:
	virtual ~AbstractPanel();
	void collapse() { this->setCollapsed(true); };
	void expand() { this->setCollapsed(false); };
	CalcSignal& clicked() { return signal_; }
	virtual void update(ChartData& chartData) { expand(); };

protected:
	AbstractPanel(std::string title, WContainerWidget* parent = 0);
	WContainerWidget* root() { return wrapper_; };
	virtual void initComponents() = 0;
	void activate(ChartData chartData) { signal_(chartData); }
	
private:
	WContainerWidget *wrapper_;
	FBspVol pdeData_;
	CalcSignal signal_;
};

#endif
