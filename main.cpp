#include <iostream>
#include <cmath>
#include <vector>
#include <set>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

using namespace std;

// Класс расчётной точки
class CalcNode {
// Класс сетки будет friend-ом точки
    friend class CalcMesh;

protected:
    // Координаты
    double x;
    double y;
    double z;
    // Выделяем шапку
    double hat;
    // Выделяем позвоночник
    double spine;
    // Угловая скорость
    double wx;
    double wy;
    double wz;

public:
    // Конструктор по умолчанию
    CalcNode() : x(0.0), y(0.0), z(0.0), hat(0.0), spine(0.0), wx(0.0), wy(0.0), wz(0.0) {
    }

    // Конструктор с указанием всех параметров
    CalcNode(double x, double y, double z, double hat, double spine, double wx)
            : x(x), y(y), z(z), hat(hat), spine(spine), wx(wx) {
    }

    // Перемещение точки
    // Вращение вокруг точки (-5, 1, 0)
    void move(double tau) {
        x += (wy * z - wz * (y - 1)) * tau;
        y += (wz * (x + 5) - wx * z) * tau;
        z += (wx * (y - 1) - wy * (x + 5)) * tau;
    }
};

// Класс элемента сетки
class Element {
// Класс сетки будет friend-ом и элемента тоже
// (и вообще будет нагло считать его просто структурой)
    friend class CalcMesh;

protected:
    // Индексы узлов, образующих этот элемент сетки
    unsigned long nodesIds[4];
};

// Класс расчётной сетки
class CalcMesh {
protected:
    // 3D-сетка из расчётных точек
    vector<CalcNode> nodes;
    vector<Element> elements;

public:
    // Конструктор сетки из заданного stl-файла
    CalcMesh(const std::vector<double> &nodesCoords, const std::vector<std::size_t> &tetrsPoints) {

        // Пройдём по узлам в модели gmsh
        nodes.resize(nodesCoords.size() / 3);
        for (unsigned int i = 0; i < nodesCoords.size() / 3; ++i) {
            // Координаты заберём из gmsh
            double pointX = nodesCoords[i * 3];
            double pointY = nodesCoords[i * 3 + 1];
            double pointZ = nodesCoords[i * 3 + 2];
            // Выделяем шапку
            double hat = pointX > 0 ? 1 : 0;
            // Выделяем позвоночник
            double spine = (pow(pointY, 2) + pow(pointZ, 2)) < 2.3 && pointX < -5 ? 1 : 0;
            // Задаём начальную угловую скорость для поворота головы
            double wx;
            if (pointX > -5) {
                wx = 1;
            } else if (pointX > -6) {
                wx = pointX + 6;
            } else {
                wx = 0;
            }
            nodes[i] = CalcNode(pointX, pointY, pointZ, hat, spine, wx);
        }

        // Пройдём по элементам в модели gmsh
        elements.resize(tetrsPoints.size() / 4);
        for (unsigned int i = 0; i < tetrsPoints.size() / 4; ++i) {
            elements[i].nodesIds[0] = tetrsPoints[i * 4] - 1;
            elements[i].nodesIds[1] = tetrsPoints[i * 4 + 1] - 1;
            elements[i].nodesIds[2] = tetrsPoints[i * 4 + 2] - 1;
            elements[i].nodesIds[3] = tetrsPoints[i * 4 + 3] - 1;
        }
    }

    // Метод задаёт угловую скорость точкам (при вращении у всех одинаковая, так что технически это честно)
    void setw(double wx, double wy, double wz) {
        for (unsigned int i = 0; i < nodes.size(); ++i) {
            // Не совсем честно использовать тут условие, функция получается не очень общая. Но мы используем.
            if (nodes[i].x > -5) {
                nodes[i].wx = wx;
                nodes[i].wy = wy;
                nodes[i].wz = wz;
            }
        }
    }

    // Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    void doTimeStep(double tau) {
        // По сути метод просто двигает все точки
        for (unsigned int i = 0; i < nodes.size(); ++i) {
            nodes[i].move(tau);
        }
    }

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    void snapshot(unsigned int snap_number) {
        // Сетка в терминах VTK
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        // Точки сетки в терминах VTK
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        // Шапка и позвоночник
        auto hat = vtkSmartPointer<vtkDoubleArray>::New();
        hat->SetName("hat");
        auto spine = vtkSmartPointer<vtkDoubleArray>::New();
        spine->SetName("spine");

        // Векторное поле на точках сетки
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        // Обходим все точки нашей расчётной сетки
        for (unsigned int i = 0; i < nodes.size(); ++i) {
            // Вставляем новую точку в сетку VTK-снапшота
            dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);

            // Добавляем значение векторного поля в этой точке
            double vx = nodes[i].wy * nodes[i].z - nodes[i].wz * (nodes[i].y - 1);
            double vy = nodes[i].wz * (nodes[i].x + 5) - nodes[i].wx * nodes[i].z;
            double vz = nodes[i].wx * (nodes[i].y - 1) - nodes[i].wy * (nodes[i].x + 5);
            double _vel[3] = {vx, vy, vz};
            vel->InsertNextTuple(_vel);

            // И значение скалярного поля тоже
            hat->InsertNextValue(nodes[i].hat);
            spine->InsertNextValue(nodes[i].spine);
        }

        // Грузим точки в сетку
        unstructuredGrid->SetPoints(dumpPoints);

        // Присоединяем векторное и скалярное поля к точкам
        unstructuredGrid->GetPointData()->AddArray(vel);
        unstructuredGrid->GetPointData()->AddArray(hat);
        unstructuredGrid->GetPointData()->AddArray(spine);

        // А теперь пишем, как наши точки объединены в тетраэдры
        for (unsigned int i = 0; i < elements.size(); ++i) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId(0, elements[i].nodesIds[0]);
            tetra->GetPointIds()->SetId(1, elements[i].nodesIds[1]);
            tetra->GetPointIds()->SetId(2, elements[i].nodesIds[2]);
            tetra->GetPointIds()->SetId(3, elements[i].nodesIds[3]);
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        // Создаём снапшот в файле с заданным именем
        string fileName = "vtk/aleut-step" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
};

int main() {
    double tau = 0.01;

    const unsigned int GMSH_TETR_CODE = 4;

    gmsh::initialize();
    gmsh::model::add("aleut");

    // Считаем STL
    try {
        gmsh::merge("fig.stl");
    } catch (...) {
        gmsh::logger::write("Could not load STL mesh: bye!");
        gmsh::finalize();
        return -1;
    }

    // Восстановим геометрию
    double angle = 20;
    bool forceParametrizablePatches = false;
    bool includeBoundary = true;
    double curveAngle = 180;
    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary, forceParametrizablePatches,
                                        curveAngle * M_PI / 180.);
    gmsh::model::mesh::createGeometry();

    // Зададим объём по считанной поверхности
    std::vector<std::pair<int, int> > s;
    gmsh::model::getEntities(s, 2);
    std::vector<int> sl;
    for (auto surf : s) sl.push_back(surf.second);
    int l = gmsh::model::geo::addSurfaceLoop(sl);
    gmsh::model::geo::addVolume({l});

    gmsh::model::geo::synchronize();

    // Зададим мелкость желаемой сетки
    int f = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(f, "F", "1");
    gmsh::model::mesh::field::setAsBackgroundMesh(f);

    // Построим сетку
    gmsh::model::mesh::generate(3);

    // Теперь извлечём из gmsh данные об узлах сетки
    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    // И данные об элементах сетки тоже извлечём, нам среди них нужны только тетраэдры, которыми залит объём
    std::vector<std::size_t> *tetrsNodesTags = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
    for (unsigned int i = 0; i < elementTypes.size(); ++i) {
        if (elementTypes[i] != GMSH_TETR_CODE)
            continue;
        tetrsNodesTags = &elementNodeTags[i];
    }

    if (tetrsNodesTags == nullptr) {
        cout << "Can not find tetra data. Exiting." << endl;
        gmsh::finalize();
        return -2;
    }

    cout << "The model has " << nodeTags.size() << " nodes and " << tetrsNodesTags->size() / 4 << " tetrs." << endl;

    // На всякий случай проверим, что номера узлов идут подряд и без пробелов
    for (int i = 0; i < nodeTags.size(); ++i) {
        // Индексация в gmsh начинается с 1, а не с нуля. Ну штош, значит так.
        assert(i == nodeTags[i] - 1);
    }
    // И ещё проверим, что в тетраэдрах что-то похожее на правду лежит.
    assert(tetrsNodesTags->size() % 4 == 0);

    CalcMesh mesh(nodesCoord, *tetrsNodesTags);

    gmsh::finalize();

    mesh.snapshot(0);

    // Глядим влево
    for (int i = 1; i <= 25; ++i) {
        mesh.doTimeStep(0.01);
        mesh.snapshot(i);
    }
    // Глядим вправо
    for (int i = 26; i <= 75; ++i) {
        mesh.doTimeStep(-0.01);
        mesh.snapshot(i);
    }

    // Голову в исходное положение
    for (int i = 76; i <= 100; ++i) {
        mesh.doTimeStep(0.01);
        mesh.snapshot(i);
    }

    mesh.setw(0, 0, 1);

    // Киваем
    for (int i = 101; i <= 110; ++i) {
        mesh.doTimeStep(0.01);
        mesh.snapshot(i);
    }
    for (int i = 111; i <= 120; ++i) {
        mesh.doTimeStep(-0.01);
        mesh.snapshot(i);
    }

    return 0;
}