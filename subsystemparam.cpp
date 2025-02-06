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

Laser SubsystemParam::getLaserData() const
{
    Laser laser;
    laser.centralWavelength = ui->lineEdit_centralWavelength->text().toDouble();
    laser.phase = ui->lineEdit_phase->text().toDouble();
    laser.pulseDuration = ui->lineEdit_pulseDuration->text().toDouble();
    return laser;
}

void SubsystemParam::on_pushButton_papram_OK_clicked()
{
    laserData = getLaserData();
    qDebug() << laserData.centralWavelength
             << laserData.phase;

}

