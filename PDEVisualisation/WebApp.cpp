#include "WebApp.h"

#include "structures.h"

#include <Wt/WBootstrapTheme>
#include <Wt/WOverlayLoadingIndicator>
#include <Wt/WText>

WebApp::WebApp(const WEnvironment& env) : WApplication(env)
{
	this->setTitle("PDE Visualisation");
	WBootstrapTheme *theme = new WBootstrapTheme(this);
	theme->setVersion(WBootstrapTheme::Version3);
	theme->setResponsive(true);
	this->setTheme(theme);
	this->useStyleSheet("resources/style.css");

	messageResourceBundle().use(appRoot() + "templates");

	wrapper_ = new WContainerWidget(root());
	root()->addStyleClass("container");

	WLoadingIndicator *loadingIndicator = new WOverlayLoadingIndicator();
	loadingIndicator->setMessage("Processing...");
	this->setLoadingIndicator(loadingIndicator);
	
	initComponents();
}

void WebApp::initComponents()
{
	wrapper_->setContentAlignment(AlignCenter);
	wrapper_->addWidget(new WText("<h1>PDE Visualisation</h1>", wrapper_));

	inputPanel_ = new InputPanel(wrapper_);
	volumePanel_ = new Chart3DPanel(wrapper_);
	slicePanel_ = new Chart2DPanel(wrapper_);

	inputPanel_->clicked().connect(boost::bind(&AbstractPanel::update, volumePanel_, _1));
	inputPanel_->clicked().connect(boost::bind(&AbstractPanel::update, slicePanel_, _1));
}


WebApp::~WebApp()
{
	my_delete(inputPanel_);
	my_delete(slicePanel_);
	my_delete(volumePanel_);
	my_delete(wrapper_);
}
