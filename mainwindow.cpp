#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_pushButton_Start_clicked()
{
    laser.centralWavelength = ui->lineEdit_centralWavelength->text().toDouble();
    laser.phase = ui->lineEdit_phase->text().toDouble();
    laser.pulseDuration = ui->lineEdit_pulseDuration->text().toDouble();
    laser.averageCountPhotons = ui->lineEdit_averageCountPh->text().toDouble();

    plotGraph(laser);
}

void MainWindow::plotGraph(const Laser &laser)
{
    // Определение параметров
    double A = 1.0;      // Амплитуда
    double omega = (3e8 / laser.centralWavelength) * 3.1415;  // Частота
    double t0 = 0.0;     // Сдвиг по времени
    double phi = 0.0;    // Фаза

    // Генерация точек графика
    QVector<double> t, E;
    // Установка диапазона для оси X, исходя из длительности импульса
    double t_min = -laser.pulseDuration;  // Начало по времени
    double t_max = laser.pulseDuration;   // Конец по времени

    for (double x = t_min; x <= t_max; x += 1e-12) {  // Диапазон X
        double y = A * cos(omega * (x - t0) + phi) * exp(-pow(x - t0, 2));
        t.append(x);
        E.append(y);
    }

    // Удаляем предыдущие графики
    ui->pulse_plot->clearGraphs();

    // Добавляем новый график
    ui->pulse_plot->addGraph();
    ui->pulse_plot->graph(0)->setData(t, E);

    // Настройки осей
    ui->pulse_plot->xAxis->setLabel("t (время)");
    ui->pulse_plot->yAxis->setLabel("E(t) (амплитуда)");
    ui->pulse_plot->xAxis->setRange(t_min, t_max);

    // Установка диапазона оси Y в зависимости от амплитуды
    double y_min = -A;  // Нижняя граница по оси Y
    double y_max = A;   // Верхняя граница по оси Y
    ui->pulse_plot->yAxis->setRange(y_min, y_max);

    // Отрисовка
    ui->pulse_plot->replot();
}

