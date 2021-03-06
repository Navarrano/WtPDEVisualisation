#ifndef ABSTRACT_PANEL_H
#define ABSTRACT_PANEL_H

#include <Wt/WPanel>
#include <Wt/WContainerWidget>
#include <string>
#include <Wt/Chart/WCartesian3DChart>
#include <Wt/WComboBox>

#include <boost/signal.hpp>

#include "structures.h"

using namespace Wt;

typedef boost::signal <void(ChartData)> CalcSignal;
typedef boost::signal <void(FBspVol, ChartData)> ExtendedCalcSignal;

class AbstractPanel : public Wt::WPanel
{
public:
	virtual ~AbstractPanel();
	void collapse() { if(!this->isCollapsed()) this->setCollapsed(true); };
	void expand() { if (this->isCollapsed()) this->setCollapsed(false); };
	CalcSignal& clicked() { return signal_; }
	ExtendedCalcSignal& extendedClicked() { return extendedSignal_; }
	virtual void update(ChartData& chartData) { expand(); };

protected:
	AbstractPanel(std::string title, WContainerWidget* parent = 0);
	WContainerWidget* root() { return wrapper_; };
	virtual void initComponents() = 0;
	void activate(ChartData& chartData) { signal_(chartData); }
	void activate(FBspVol& vol, ChartData& chartData) { extendedSignal_(vol, chartData); }
private:
	WContainerWidget *wrapper_;
	FBspVol pdeData_;
	CalcSignal signal_;
	ExtendedCalcSignal extendedSignal_;
};

#endif
