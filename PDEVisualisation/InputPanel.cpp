#include "InputPanel.h"

#include "global_methods.h"

#include <Wt/WText>
#include <Wt/WGroupBox>
#include <Wt/WLabel>
#include <Wt/WTemplate>

InputPanel::InputPanel(WContainerWidget* parent) : AbstractPanel("Input data", parent)
{
	initComponents();
	initEvents();

	//std::ifstream ifs("F:/LaplaceSoln2.dat");
	//ifs >> pdeData_;
}

void InputPanel::initComponents()
{
	WGroupBox* pdeBox = new WGroupBox("Partial derivative equation", root());
	pdeBox->addStyleClass("main-legend");

	ordersBox_ = new OrdersGroupBox(pdeBox);
	limitsBox_ = new LimitsGroupBox(pdeBox);
	boundariesBox_ = new BoundariesGroupBox(pdeBox);
	rhsGroupBox_ = new RhsGroupBox(pdeBox);
	settingsBox_ = new SystemSettingsBox(root());

	new WText("<hr/>", root());
	calculateButton_ = new WPushButton("Calculate", root());
}

void InputPanel::initEvents()
{
	calculateButton_->clicked().connect(std::bind([=](){
		ChartData chartData;
		unsigned density = settingsBox_->density();
		pdeData_ = calculatePDE(ordersBox_->getOrders(), limitsBox_->getLimits(), boundariesBox_->boundaries(), rhsGroupBox_->rhs());
		chartData.xStart = pdeData_.GetLeftLimitU();
		chartData.xEnd = pdeData_.GetRightLimitU();
		chartData.yStart = pdeData_.GetLeftLimitV();
		chartData.yEnd = pdeData_.GetRightLimitV();
		chartData.zStart = pdeData_.GetLeftLimitW();
		chartData.zEnd = pdeData_.GetRightLimitW();
		chartData.xSize = (chartData.xEnd - chartData.xStart) * density + 1;
		chartData.ySize = (chartData.yEnd - chartData.yStart) * density + 1;
		chartData.zSize = (chartData.zEnd - chartData.zStart) * density + 1;
		chartData.totalNbPts = chartData.xSize * chartData.ySize * chartData.zSize;
		chartData.points = pdeData_.ComputePoints(chartData.xSize - 1, chartData.ySize - 1, chartData.zSize - 1);
		chartData.step = settingsBox_->stepSize();
		activate(chartData);
	}));

	limitsBox_->valid().connect(std::bind([=](){
		if (!boundariesBox_->enabled())
		{
			boundariesBox_->enableAll();
		}
		boundariesBox_->update(limitsBox_->getLimits());
	}));
	limitsBox_->invalid().connect(std::bind([=](){
		if (boundariesBox_->enabled())
		{
			boundariesBox_->disableAll();
			rhsGroupBox_->disable();
			calculateButton_->setDisabled(true);
		}
	}));
	boundariesBox_->valid().connect(std::bind([=](){
		rhsGroupBox_->enable();
	}));
	boundariesBox_->invalid().connect(std::bind([=](){
		rhsGroupBox_->disable();
		calculateButton_->setDisabled(true);
	}));
	rhsGroupBox_->valid().connect(std::bind([=](){
		calculateButton_->setDisabled(false);
	}));
	rhsGroupBox_->invalid().connect(std::bind([=](){
		calculateButton_->setDisabled(true);
	}));
}

FBspVol InputPanel::preperVolData(SurfaceOrders orders, Limits limits)
{
	short orderx = orders.xOrder;
	short ordery = orders.yOrder;
	short orderz = orders.zOrder;
	short segments = 5;

	double stepX = (limits.xRightLimit - limits.xLeftLimit) / (double)segments;
	double stepY = (limits.yRightLimit - limits.yLeftLimit) / (double)segments;
	double stepZ = (limits.zRightLimit - limits.zLeftLimit) / (double)segments;

	short numx = orderx + segments - 1;
	short numy = ordery + segments - 1;
	short numz = orderz + segments - 1;

	// knot sets
	Vector<double> knotX(numx + orderx), knotY(numy + ordery), knotZ(numz + orderz);

	for (int i = 0; i<orderx; i++) knotX[i] = limits.xLeftLimit;
	for (int i = orderx; i<numx; i++) knotX[i] = limits.xLeftLimit + (i - orderx + 1)*stepX;
	for (int i = numx; i<numx + orderx; i++) knotX[i] = limits.xRightLimit;

	for (int i = 0; i<ordery; i++) knotY[i] = limits.yLeftLimit;
	for (int i = ordery; i<numy; i++) knotY[i] = limits.yLeftLimit + (i - ordery + 1)*stepY;
	for (int i = numy; i<numy + ordery; i++) knotY[i] = limits.yRightLimit;

	for (int i = 0; i<orderz; i++) knotZ[i] = limits.zLeftLimit;
	for (int i = orderz; i<numz; i++) knotZ[i] = limits.zLeftLimit + (i - orderz + 1)*stepZ;
	for (int i = numz; i<numz + orderz; i++) knotZ[i] = limits.zRightLimit;

	Matrix3D<double> cpts(numx, numy, numz, 0.0);

	return FBspVol(cpts, knotX, knotY, knotZ, orderx, ordery, orderz, numx, numy, numz);
}

FBspVol InputPanel::calculatePDE(SurfaceOrders surfaceOrders, Limits limits, std::array<PolynomialData, 6> boundaries, PolynomialData rhs)
{
	FBspVol b = preperVolData(surfaceOrders, limits);

	int num_pointforce = 0, num_distforce = 1, num_supports = 0;

	std::vector<int> nums(40);

	nums[0] = num_pointforce;
	nums[1] = num_distforce;
	nums[2] = num_supports;

	std::vector<double> pointforce(4 * num_pointforce);
	std::vector<double> distforce(6 * num_distforce);
	std::vector<double> supports(4 * num_supports);
	std::vector<int> orders(3 * num_distforce);

	for (int i = 0; i < 4 * num_pointforce; i++) pointforce[i] = 0;
	distforce[0] = limits.xLeftLimit;
	distforce[1] = limits.xRightLimit;
	distforce[2] = limits.yLeftLimit;
	distforce[3] = limits.yRightLimit;
	distforce[4] = limits.zLeftLimit;
	distforce[5] = limits.zRightLimit;

	for (int i = 0; i < 3 * num_distforce; i++) orders[i] = rhs.order_;

	Vector<Matrix3D<double> > vol(num_distforce);
	int ind = 0;

	if (num_distforce > 0) {
		for (int i = 0; i<num_distforce; i++) {
			Matrix3D<double> voldata(orders[ind], orders[ind + 1], orders[ind + 2]);
			int l = 0;
			for (int k = 0; k<orders[ind + 2]; k++)
				for (int i = 0; i<orders[ind]; i++)
					for (int j = 0; j<orders[ind + 1]; j++) voldata[k][i][j] = rhs.coeffs_[l++];
			ind = ind + 3;
			vol[i] = voldata;
		}
	}

	for (int i = 0; i<4 * num_supports; i++) supports[i] = 0;
	double L = limits.xRightLimit - limits.xLeftLimit;
	double W = limits.yRightLimit - limits.yLeftLimit;
	double H = limits.zRightLimit - limits.zLeftLimit;

	// boundary conditions

	Vector<int> natgeom(6);
	for (int i = 0; i < 6; i++) natgeom[i] = 0;

	fillVector(nums, boundaries[2].order_, 3, 2);
	fillVector(nums, 0, 5, 4);
	fillVector(nums, boundaries[3].order_, 9, 2);
	fillVector(nums, 0, 11, 4);
	fillVector(nums, boundaries[5].order_, 15, 2);
	fillVector(nums, 0, 17, 4);
	fillVector(nums, boundaries[4].order_, 21, 2);
	fillVector(nums, 0, 23, 4);
	fillVector(nums, boundaries[1].order_, 27, 2);
	fillVector(nums, 0, 29, 4);
	fillVector(nums, boundaries[0].order_, 33, 2);
	fillVector(nums, 0, 35, 4);

	int num = Math::max1(Math::max1(b.GetNumU(), b.GetNumV()), b.GetNumW());
	Matrix3D<double> bound(num, num, 18, 0.0);

	
	// boundary 1
	int k = 0;
	for (int i = 0; i<nums[3]; i++)
		for (int j = 0; j<nums[4]; j++) {
			bound[0][i][j] = boundaries[2].coeffs_[k++];
		}
	for (int i = 0; i<nums[5]; i++)
		for (int j = 0; j<nums[6]; j++) {
			bound[1][i][j] = 0;
		}
	for (int i = 0; i<nums[7]; i++)
		for (int j = 0; j<nums[8]; j++) {
			bound[2][i][j] = 0;
		}

	// boundary 2
	k = 0;
	for (int i = 0; i<nums[9]; i++)
		for (int j = 0; j<nums[10]; j++) {
			bound[3][i][j] = boundaries[3].coeffs_[k++];
		}
	for (int i = 0; i<nums[11]; i++)
		for (int j = 0; j<nums[12]; j++) {
			bound[4][i][j] = 0;
		}
	for (int i = 0; i<nums[13]; i++)
		for (int j = 0; j<nums[14]; j++) {
			bound[5][i][j] = 0;
		}

	// boundary 3
	k = 0;
	for (int i = 0; i<nums[15]; i++)
		for (int j = 0; j<nums[16]; j++) {
			bound[6][i][j] = boundaries[5].coeffs_[k++];
		}
	for (int i = 0; i<nums[17]; i++)
		for (int j = 0; j<nums[18]; j++) {
			bound[7][i][j] = 0;
		}
	for (int i = 0; i<nums[19]; i++)
		for (int j = 0; j<nums[20]; j++) {
			bound[8][i][j] = 0;
		}

	// boundary 4
	k = 0;
	for (int i = 0; i<nums[21]; i++)
		for (int j = 0; j<nums[22]; j++) {
			bound[9][i][j] = boundaries[4].coeffs_[k++];
		}
	for (int i = 0; i<nums[23]; i++)
		for (int j = 0; j<nums[24]; j++) {
			bound[10][i][j] = 0;
		}
	for (int i = 0; i<nums[25]; i++)
		for (int j = 0; j<nums[26]; j++) {
			bound[11][i][j] = 0;
		}

	// boundary 5
	k = 0;
	for (int i = 0; i<nums[27]; i++)
		for (int j = 0; j<nums[28]; j++) {
			bound[12][i][j] = boundaries[1].coeffs_[k++];
		}
	for (int i = 0; i<nums[29]; i++)
		for (int j = 0; j<nums[30]; j++) {
			bound[13][i][j] = 0;
		}
	for (int i = 0; i<nums[31]; i++)
		for (int j = 0; j<nums[32]; j++) {
			bound[14][i][j] = 0;
		}

	// boundary 6
	k = 0;
	for (int i = 0; i<nums[33]; i++)
		for (int j = 0; j<nums[34]; j++) {
			bound[15][i][j] = boundaries[0].coeffs_[k++];
		}
	for (int i = 0; i<nums[35]; i++)
		for (int j = 0; j<nums[36]; j++) {
			bound[16][i][j] = 0;
		}
	for (int i = 0; i<nums[37]; i++)
		for (int j = 0; j<nums[38]; j++) {
			bound[17][i][j] = 0;
		}
	nums[39] = 1;

	if (rhs.order_ == 1)
		return b.ComputeSolution(L, W, H, Vector<int>(nums), Vector<double>(supports), Vector<double>(pointforce), Vector<double>(distforce), Vector<int>(orders), vol, natgeom, bound);
	else
		return b.ComputeSolution1(L, W, H, Vector<int>(nums), Vector<double>(supports), Vector<double>(pointforce), Vector<double>(distforce), Vector<int>(orders), vol, natgeom, bound);
	return b;
}

void InputPanel::fillVector(std::vector<int>& nums, int val, size_t start, size_t lenght)
{
	for (int i = start; i < start + lenght; i++)
	{
		nums[i] = val;
	}
}

InputPanel::~InputPanel()
{
	my_delete(ordersBox_);
	my_delete(limitsBox_);
	my_delete(boundariesBox_);
	my_delete(rhsGroupBox_);
	my_delete(settingsBox_);

	my_delete(calculateButton_);
}