#include "AbstractPanel.h"

#include <Wt/WAnimation>

AbstractPanel::AbstractPanel(std::string title, WContainerWidget* parent) : WPanel(parent)
{
	this->setTitle(title);
	this->setCollapsible(true);
	this->addStyleClass("panel-primary");
	this->setMargin(30, Bottom);

	WAnimation animation(WAnimation::SlideInFromTop | WAnimation::Fade, WAnimation::EaseOut, 200);
	this->setAnimation(animation);

	wrapper_ = new WContainerWidget();
	this->setCentralWidget(wrapper_);
}

AbstractPanel::~AbstractPanel()
{
	delete wrapper_;
}
