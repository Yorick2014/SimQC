#ifndef SUBSYSTEMPARAM_H
#define SUBSYSTEMPARAM_H

#include <QWidget>

namespace Ui {
class SubsystemParam;
}

class SubsystemParam : public QWidget
{
    Q_OBJECT

public:
    explicit SubsystemParam(QWidget *parent = nullptr);
    ~SubsystemParam();

private:
    Ui::SubsystemParam *ui;
};

#endif // SUBSYSTEMPARAM_H
