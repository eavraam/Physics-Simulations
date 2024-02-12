#ifndef WIDGETFLUID_H
#define WIDGETFLUID_H

#include <QWidget>

namespace Ui {
class WidgetFluid;
}

class WidgetFluid : public QWidget
{
    Q_OBJECT
public:
    explicit WidgetFluid(QWidget *parent = nullptr);
    ~WidgetFluid();

    double getGravity()    const;

signals:
    void updatedParameters();

private:
    Ui::WidgetFluid *ui;
};

#endif // WIDGETFLUID_H
