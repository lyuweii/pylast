#pragma once
#include "ArrayEvent.hh"
#include "EventSource.hh"
#include <arrow/api.h>
#include <memory>
#include <optional>
#include <string>
#include <thread>
#include <vector>

namespace df {
    struct SimulationData {
        SimulationData(size_t reserve_rows = 0) { init_builders(reserve_rows); }

        static std::shared_ptr<arrow::Schema> make_schema() {
            return arrow::schema(
                {arrow::field("event_id", arrow::int32()),
                 arrow::field("true_energy", arrow::float64()),
                 arrow::field("true_alt", arrow::float64()),
                 arrow::field("true_az", arrow::float64()),
                 arrow::field("true_core_x", arrow::float64()),
                 arrow::field("true_core_y", arrow::float64()),
                 arrow::field("h_first_int", arrow::float64()),
                 arrow::field("x_max", arrow::float64()),
                 arrow::field("h_max", arrow::float64()),
                 arrow::field("starting_grammage", arrow::float64())});
        }

        void write_event(const ArrayEvent &event) {
            const auto &shower = event.simulation->shower;
            event_id.Append(event.event_id);
            true_energy.Append(shower.energy);
            true_alt.Append(shower.alt);
            true_az.Append(shower.az);
            true_core_x.Append(shower.core_x);
            true_core_y.Append(shower.core_y);
            h_first_int.Append(shower.h_first_int);
            x_max.Append(shower.x_max);
            h_max.Append(shower.h_max);
            starting_grammage.Append(shower.starting_grammage);
        }

        void init_builders(size_t reserve_rows) {
            if (reserve_rows == 0)
                reserve_rows = 1000; // 默认预留1000行
            event_id.Reserve(reserve_rows);
            true_energy.Reserve(reserve_rows);
            true_alt.Reserve(reserve_rows);
            true_az.Reserve(reserve_rows);
            true_core_x.Reserve(reserve_rows);
            true_core_y.Reserve(reserve_rows);
            h_first_int.Reserve(reserve_rows);
            x_max.Reserve(reserve_rows);
            h_max.Reserve(reserve_rows);
            starting_grammage.Reserve(reserve_rows);
        }

        void
        finish_builders(std::vector<std::shared_ptr<arrow::Array>> &array) {
            arrow::ArrayVector arrs(10);
            event_id.Finish(&arrs[0]);
            true_energy.Finish(&arrs[1]);
            true_alt.Finish(&arrs[2]);
            true_az.Finish(&arrs[3]);
            true_core_x.Finish(&arrs[4]);
            true_core_y.Finish(&arrs[5]);
            h_first_int.Finish(&arrs[6]);
            x_max.Finish(&arrs[7]);
            h_max.Finish(&arrs[8]);
            starting_grammage.Finish(&arrs[9]);
            array = std::move(arrs);
        }

        arrow::Int32Builder event_id;
        arrow::DoubleBuilder true_energy, true_alt, true_az;
        arrow::DoubleBuilder true_core_x, true_core_y;
        arrow::DoubleBuilder h_first_int, x_max, h_max, starting_grammage;
    };

    // ------------------------ ReconstructorData ------------------------
    struct ReconstructorData {
        ReconstructorData(size_t reserve_rows = 0) {
            init_builders(reserve_rows);
        }

        static std::shared_ptr<arrow::Schema> make_schema() {
            return arrow::schema(
                {arrow::field("reconstructor_name", arrow::utf8()),
                 arrow::field("event_id", arrow::int32()),
                 arrow::field("is_valid", arrow::boolean()),
                 arrow::field("rec_alt", arrow::float64()),
                 arrow::field("rec_alt_uncertainty", arrow::float64()),
                 arrow::field("rec_az", arrow::float64()),
                 arrow::field("rec_az_uncertainty", arrow::float64()),
                 arrow::field("rec_core_x", arrow::float64()),
                 arrow::field("rec_core_y", arrow::float64()),
                 arrow::field("rec_core_position_error", arrow::float64()),
                 arrow::field("tilted_core_x", arrow::float64()),
                 arrow::field("tilted_core_y", arrow::float64()),
                 arrow::field("tilted_core_uncertainty_x", arrow::float64()),
                 arrow::field("tilted_core_uncertainty_y", arrow::float64()),
                 arrow::field("hmax", arrow::float64()),
                 arrow::field("direction_error", arrow::float64())});
        }

        void write_event(const ArrayEvent &event) {
            const auto &geometry = event.dl2->geometry;
            for (const auto &[name, recon] : geometry) {
                std::cout << "recon.alt:" << recon.alt<<std::endl;
                reconstructor_name.Append(name);
                event_id.Append(event.event_id);
                is_valid.Append(recon.is_valid);
                rec_alt.Append(recon.alt);
                rec_alt_uncertainty.Append(recon.alt_uncertainty);
                rec_az.Append(recon.az);
                rec_az_uncertainty.Append(recon.az_uncertainty);
                rec_core_x.Append(recon.core_x);
                rec_core_y.Append(recon.core_y);
                rec_core_position_error.Append(recon.core_pos_error);
                tilted_core_x.Append(recon.tilted_core_x);
                tilted_core_y.Append(recon.tilted_core_y);
                tilted_core_uncertainty_x.Append(
                    recon.tilted_core_uncertainty_x);
                tilted_core_uncertainty_y.Append(
                    recon.tilted_core_uncertainty_y);
                hmax.Append(recon.hmax);
                direction_error.Append(recon.direction_error);
            }
        }

        void init_builders(size_t reserve_rows) {
            if (reserve_rows == 0)
                reserve_rows = 1000; // 默认预留1000行
            reconstructor_name.Reserve(reserve_rows);
            event_id.Reserve(reserve_rows);
            is_valid.Reserve(reserve_rows);
            rec_alt.Reserve(reserve_rows);
            rec_alt_uncertainty.Reserve(reserve_rows);
            rec_az.Reserve(reserve_rows);
            rec_az_uncertainty.Reserve(reserve_rows);
            rec_core_x.Reserve(reserve_rows);
            rec_core_y.Reserve(reserve_rows);
            rec_core_position_error.Reserve(reserve_rows);
            tilted_core_x.Reserve(reserve_rows);
            tilted_core_y.Reserve(reserve_rows);
            tilted_core_uncertainty_x.Reserve(reserve_rows);
            tilted_core_uncertainty_y.Reserve(reserve_rows);
            hmax.Reserve(reserve_rows);
            direction_error.Reserve(reserve_rows);
        }

        void
        finish_builders(std::vector<std::shared_ptr<arrow::Array>> &array) {
            arrow::ArrayVector arrs(16);
            reconstructor_name.Finish(&arrs[0]);
            event_id.Finish(&arrs[1]);
            is_valid.Finish(&arrs[2]);
            rec_alt.Finish(&arrs[3]);
            rec_alt_uncertainty.Finish(&arrs[4]);
            rec_az.Finish(&arrs[5]);
            rec_az_uncertainty.Finish(&arrs[6]);
            rec_core_x.Finish(&arrs[7]);
            rec_core_y.Finish(&arrs[8]);
            rec_core_position_error.Finish(&arrs[9]);
            tilted_core_x.Finish(&arrs[10]);
            tilted_core_y.Finish(&arrs[11]);
            tilted_core_uncertainty_x.Finish(&arrs[12]);
            tilted_core_uncertainty_y.Finish(&arrs[13]);
            hmax.Finish(&arrs[14]);
            direction_error.Finish(&arrs[15]);
            array = std::move(arrs);
        }

        arrow::StringBuilder reconstructor_name;
        arrow::Int32Builder event_id;
        arrow::BooleanBuilder is_valid;
        arrow::DoubleBuilder rec_alt, rec_alt_uncertainty;
        arrow::DoubleBuilder rec_az, rec_az_uncertainty;
        arrow::DoubleBuilder rec_core_x, rec_core_y;
        arrow::DoubleBuilder rec_core_position_error;
        arrow::DoubleBuilder tilted_core_x, tilted_core_y;
        arrow::DoubleBuilder tilted_core_uncertainty_x,
            tilted_core_uncertainty_y;
        arrow::DoubleBuilder hmax, direction_error;
    };

    struct TelescopeData {
        TelescopeData(size_t reserve_rows = 0) { init_builders(reserve_rows); }

        static std::shared_ptr<arrow::Schema> make_schema() {
            return arrow::schema({
                arrow::field("event_id", arrow::int32()),
                arrow::field("tel_id", arrow::int32()),
                arrow::field("hillas_length", arrow::float64()),
                arrow::field("hillas_width", arrow::float64()),
                arrow::field("hillas_psi", arrow::float64()),
                arrow::field("hillas_x", arrow::float64()),
                arrow::field("hillas_y", arrow::float64()),
                arrow::field("hillas_intensity", arrow::float64()),
            });
        }

        void write_event(const ArrayEvent &event) {
            if (!event.dl1.has_value()){
                return;
            }
            const auto &dl1 = event.dl1;            
            std::cout << dl1->tels.size()<<"\n";
            for (const auto &[index, tel] : dl1->tels) {
                event_id.Append(event.event_id);
                tel_id.Append(index);
                hillas_intensity.Append(tel->image_parameters.hillas.intensity);
                hillas_length.Append(tel->image_parameters.hillas.length);
                hillas_width.Append((tel->image_parameters.hillas.width));
                hillas_psi.Append(tel->image_parameters.hillas.psi);
                hillas_x.Append(tel->image_parameters.hillas.x);
                hillas_y.Append(tel->image_parameters.hillas.y);
            }
        }

        void init_builders(size_t reserve_rows) {
            if (reserve_rows == 0)
                reserve_rows = 1000; // 默认预留1000行
            event_id.Reserve(reserve_rows);
            tel_id.Reserve(reserve_rows);
            hillas_intensity.Reserve(reserve_rows);
            hillas_length.Reserve(reserve_rows);
            hillas_width.Reserve(reserve_rows);
            hillas_x.Reserve(reserve_rows);
            hillas_y.Reserve(reserve_rows);
            hillas_psi.Reserve(reserve_rows);
        }

        void
        finish_builders(std::vector<std::shared_ptr<arrow::Array>> &array) {
            arrow::ArrayVector arrs(8);
            event_id.Finish(&arrs[0]);
            tel_id.Finish(&arrs[1]);
            hillas_length.Finish(&arrs[2]);
            hillas_width.Finish(&arrs[3]);
            hillas_psi.Finish(&arrs[4]);
            hillas_x.Finish(&arrs[5]);
            hillas_y.Finish(&arrs[6]);
            hillas_intensity.Finish(&arrs[7]);
            array = std::move(arrs);
        }

        arrow::Int32Builder event_id, tel_id;
        arrow::DoubleBuilder hillas_length, hillas_width, hillas_psi;
        arrow::DoubleBuilder hillas_x, hillas_y;
        arrow::DoubleBuilder hillas_intensity;
    };

    // ------------------------ TableMaker ------------------------
    template <typename TableType>
    class TableMaker {
    public:
        TableMaker(size_t reserve_rows = 0) : _data(reserve_rows) {
            _schema = TableType::make_schema();
            _arrays.resize(_schema->num_fields());
        }

        void write_event(const ArrayEvent &event) { _data.write_event(event); }

        std::shared_ptr<arrow::Table> make_table() {
            _data.finish_builders(_arrays);
            return arrow::Table::Make(_schema, _arrays);
        }

    private:
        TableType _data;
        std::vector<std::shared_ptr<arrow::Array>> _arrays;
        std::shared_ptr<arrow::Schema> _schema;
    };

    struct DataTable {
        std::shared_ptr<arrow::Table> simulation_table;
        std::shared_ptr<arrow::Table> reconstructor_table;
        std::shared_ptr<arrow::Table> telescope_table;
    };

    class DataFrameMaker {
    public:
        DataFrameMaker(EventSource &source, size_t reserve_rows = 100000)
            : _source(source), _simulation_table_maker(reserve_rows),
              _reconstructor_table_maker(reserve_rows),
              _telscope_table_maker(reserve_rows){};

        // std::shared_ptr<arrow::Table> operator()();
        DataTable operator()();

        auto build_table() {
            for (const auto &event : _source) {
                _simulation_table_maker.write_event(event);
                _reconstructor_table_maker.write_event(event);
                _telscope_table_maker.write_event(event);
            }
            auto simulation_table = _simulation_table_maker.make_table();
            auto reconstructor_table = _reconstructor_table_maker.make_table();
            auto telescope_table = _telscope_table_maker.make_table();
            return DataTable{simulation_table, reconstructor_table,
                             telescope_table};
        }

    private:
        EventSource &_source;
        TableMaker<SimulationData> _simulation_table_maker;
        TableMaker<ReconstructorData> _reconstructor_table_maker;
        TableMaker<TelescopeData> _telscope_table_maker;
    };

    class MultiThreadedDataFrameMaker {
    public:
        MultiThreadedDataFrameMaker(EventSource &source, size_t num_threads = 4)
            : _source(source), _num_threads(num_threads) {}

        std::tuple<std::shared_ptr<arrow::Table>, std::shared_ptr<arrow::Table>>
        build_table() {
            // 拆分任务
            size_t total_events = _source.max_events;
            size_t batch_size =
                (total_events + _num_threads - 1) / _num_threads;

            std::vector<std::thread> threads;
            std::vector<std::vector<std::shared_ptr<arrow::Array>>>
                simulation_arrays(_num_threads);
            std::vector<std::vector<std::shared_ptr<arrow::Array>>>
                reconstructor_arrays(_num_threads);

            // 启动线程
            for (size_t t = 0; t < _num_threads; ++t) {
                threads.emplace_back([this, t, batch_size, &simulation_arrays,
                                      &reconstructor_arrays]() {
                    size_t start_idx = t * batch_size;
                    size_t max_events = _source.max_events;
                    size_t end_idx =
                        std::min(start_idx + batch_size, max_events);

                    // 每个线程单独 Builder
                    SimulationData sim_data(end_idx - start_idx);
                    ReconstructorData rec_data(end_idx - start_idx);

                    for (size_t i = start_idx; i < end_idx; ++i) {
                        sim_data.write_event(_source[i]);
                        rec_data.write_event(_source[i]);
                    }

                    sim_data.finish_builders(simulation_arrays[t]);
                    rec_data.finish_builders(reconstructor_arrays[t]);
                });
            }

            // 等待线程完成
            for (auto &th : threads)
                th.join();

            // 合并每个线程生成的 Array
            std::vector<std::shared_ptr<arrow::Array>> merged_sim_arrays;
            std::vector<std::shared_ptr<arrow::Array>> merged_rec_arrays;

            merged_sim_arrays.resize(simulation_arrays[0].size());
            merged_rec_arrays.resize(reconstructor_arrays[0].size());

            for (size_t col = 0; col < merged_sim_arrays.size(); ++col) {
                merged_sim_arrays[col] =
                    *arrow::Concatenate(simulation_arrays[col]);
                merged_rec_arrays[col] =
                    *arrow::Concatenate(reconstructor_arrays[col]);
            }

            auto simulation_table = arrow::Table::Make(
                SimulationData::make_schema(), merged_sim_arrays);
            auto reconstructor_table = arrow::Table::Make(
                ReconstructorData::make_schema(), merged_rec_arrays);

            return {simulation_table, reconstructor_table};
        }

    private:
        EventSource &_source;
        size_t _num_threads;
    };

} // namespace df