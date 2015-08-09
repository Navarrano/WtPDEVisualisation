#ifndef CHART_3D_PANEL_H
#define CHART_3D_PANEL_H

#include <Wt/WGroupBox>
#include <Wt/WComboBox>
#include <Wt/WRadioButton>
#include <Wt/WButtonGroup>
#include <Wt/WSlider>
#include <Wt/WCheckBox>

#include <boost/signals.hpp>

#include "BaseChart.h"
#include "AbstractPanel.h"
#include "Models.h"

class SettingsBox;
class DisplayRadioButtons;
class ClippingGroupBox;

struct VolumeSettings 
{
	DisplayRate xRate, yRate, zRate;
	short generalSize, skippedSize;
	double xStart, xEnd, yStart, yEnd, zStart, zEnd;
};

class Chart3DPanel :
	public AbstractPanel
{
public:
	Chart3DPanel(WContainerWidget *parent = 0);
	virtual void update(ChartData& chartData);
	virtual ~Chart3DPanel();
private:
	void initComponents();
	void initChartLook();
	void initEvents();

	BaseChart* volumeChart_;
	SettingsBox* settingsBox_;
	VolumeData* volumeData_;
	WScatterData* scatter_;
};

////// SETTINGS BOX ////////

class SettingsBox : public WContainerWidget
{
public:
	SettingsBox(WContainerWidget *parent = 0);
	void update(ChartData& chartData);
	boost::signal<void(DisplayRate, DisplayRate, DisplayRate)>& changedRate() { return rateChangeSignal_; };
	boost::signal<void(VolumeSettings)>& changed() { return settingsChangedSignal_; };
	boost::signal<void(Limits)>& clip();
	virtual ~SettingsBox();
	WCheckBox* getCheckBox() { return showColorMap_; };
private:
	void initComponenets();
	void initEvents();
	void activateControls();
	void initComboBox(WComboBox* cb, int index);
	VolumeSettings generateSettings();

	DisplayRadioButtons* xRate_;
	DisplayRadioButtons* yRate_;
	DisplayRadioButtons* zRate_;
	ClippingGroupBox* clippingBox_;
	WComboBox *general_;
	WComboBox *skipped_;
	WCheckBox *showColorMap_;
	boost::signal<void(DisplayRate, DisplayRate, DisplayRate)> rateChangeSignal_;
	boost::signal<void(VolumeSettings)> settingsChangedSignal_;
};

//////// DISPLAY RADIO BUTTONS ////////

class DisplayRadioButtons : public WContainerWidget
{
public:
	DisplayRadioButtons(char* label, WContainerWidget *parent = 0);
	boost::signal<void()>& changed() { return radioChangedSignal_; };
	DisplayRate getRate();
	void activate();
	virtual ~DisplayRadioButtons();
private:
	void initComponents(char* label);
	void initEvents();

	WButtonGroup *group_;
	WRadioButton *everyPoint_;
	WRadioButton *everySecondPoint_;
	WRadioButton *everyFifthPoint_;
	WRadioButton *everyTenthPoint_;
	boost::signal<void()> radioChangedSignal_;
};

//////// CLIPPING SLIDERS ////////

class ClippingAxisSliders : public WContainerWidget
{
public:
	ClippingAxisSliders(char* label, WContainerWidget* parent = 0);
	void update(double min, double step, unsigned nbPts);
	virtual ~ClippingAxisSliders();
	double getLeftLimit();
	double getRightLimit();
private:
	void initComponents();
	void initEvents();
	WString tickToTextValue(int tick);
	WSlider* initSlider();

	WSlider* leftLimitSlider_;
	WSlider* rightLimitSlider_;
	WText* leftLimitText_;
	WText* rightLimitText_;

	double min_;
	double step_;
};


//////// CLIPPING GROUP BOX ///////////

class ClippingGroupBox : public WGroupBox
{
public:
	ClippingGroupBox(WContainerWidget* parent = 0);
	void update(ChartData& chartData);
	boost::signal<void(Limits)>& clicked() { return clipButtonSignal_; };
	virtual ~ClippingGroupBox();
private:
	void initComponents();
	void initEvents();

	ClippingAxisSliders* uSliders_;
	ClippingAxisSliders* vSliders_;
	ClippingAxisSliders* wSliders_;
	WPushButton* clip_;
	boost::signal<void(Limits)> clipButtonSignal_;
};

#endif
