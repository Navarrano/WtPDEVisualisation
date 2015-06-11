#include "ViewButtonsWidget.h"


ViewButtonsWidget::ViewButtonsWidget(WContainerWidget* parent) : WContainerWidget(parent)
{
	perspectiveView_ = new WPushButton("Perspective View", this);
	topView_ = new WPushButton("Top View", this);
	frontView_ = new WPushButton("Front view", this);
	sideView_ = new WPushButton("Side View", this);

	perspectiveView_->clicked().connect(std::bind([=]() {
		viewChangeSignal_(ViewType::PERSPECTIVE);
	}));
	topView_->clicked().connect(std::bind([=]() {
		viewChangeSignal_(ViewType::TOP);
	}));
	sideView_->clicked().connect(std::bind([=]() {
		viewChangeSignal_(ViewType::SIDE);
	}));
	frontView_->clicked().connect(std::bind([=]() {
		viewChangeSignal_(ViewType::FRONT);
	}));
}

ViewButtonsWidget::~ViewButtonsWidget()
{
	delete perspectiveView_;
	delete topView_;
	delete sideView_;
	delete frontView_;
}
