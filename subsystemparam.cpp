#include "subsystemparam.h"
#include "ui_subsystemparam.h"

SubsystemParam::SubsystemParam(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::SubsystemParam)
{
    ui->setupUi(this);
}

SubsystemParam::~SubsystemParam()
{
    delete ui;
}
