#include "widgetfluid.h"
#include "ui_widgetfluid.h"

WidgetFluid::WidgetFluid(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::WidgetFluid)
{
    ui->setupUi(this);

    connect(ui->btnUpdate, &QPushButton::clicked, this,
            [=] (void) { emit updatedParameters(); });
}

WidgetFluid::~WidgetFluid()
{
    delete ui;
}

double WidgetFluid::getGravity() const {
    return ui->gravity->value();
}
