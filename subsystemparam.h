#ifndef SUBSYSTEMPARAM_H
#define SUBSYSTEMPARAM_H

#include <QWidget>
#include "components.h"

namespace Ui {
class SubsystemParam;
}

class SubsystemParam : public QWidget
{
    Q_OBJECT

public:
    explicit SubsystemParam(QWidget *parent = nullptr);
    ~SubsystemParam();

    Laser getLaserData() const;

private slots:
    void on_pushButton_papram_OK_clicked();

private:
    Ui::SubsystemParam *ui;
    Laser laserData;
};

#endif // SUBSYSTEMPARAM_H
