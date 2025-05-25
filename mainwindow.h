#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "components.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:

    void on_pushButton_Start_clicked();

private:
    Ui::MainWindow *ui;
    Laser laser;
    QuantumChannel quantumChannel;

    double generate_random_0_to_1();
    void plotGraph(const Laser &laser);
    void plotTimeDomain(const Laser &laser);
    void plotGenKeys(const Laser &laser);
};
#endif // MAINWINDOW_H
