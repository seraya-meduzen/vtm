#include <mavsdk/mavsdk.h>
#include <mavsdk/plugins/action/action.h>
#include <mavsdk/plugins/mission/mission.h>
#include <mavsdk/plugins/telemetry/telemetry.h>
#include <mavsdk/plugins/calibration/calibration.h>

#include <chrono>
#include <functional>
#include <future>
#include <iostream>
#include <thread>
#include <fstream>
#include <cmath>

using namespace mavsdk;
using std::chrono::seconds;
using std::this_thread::sleep_for;
using std::isnan;

static std::function<void(Calibration::Result, Calibration::ProgressData)>
create_calibration_callback(std::promise<void>&);

static void calibrate_accelerometer(Calibration&);
static void calibrate_gyro(Calibration&);
static void calibrate_magnetometer(Calibration&);



Mission::MissionItem make_mission_item(
    double latitude_deg,
    double longitude_deg,
    float relative_altitude_m,
    float speed_m_s,
    bool is_fly_through,
    float gimbal_pitch_deg,
    float gimbal_yaw_deg,
    Mission::MissionItem::CameraAction camera_action)
{
    Mission::MissionItem new_item{};
    new_item.latitude_deg = latitude_deg;
    new_item.longitude_deg = longitude_deg;
    new_item.relative_altitude_m = relative_altitude_m;
    new_item.speed_m_s = speed_m_s;
    new_item.is_fly_through = is_fly_through;
    new_item.gimbal_pitch_deg = gimbal_pitch_deg;
    new_item.gimbal_yaw_deg = gimbal_yaw_deg;
    new_item.camera_action = camera_action;
    return new_item;
}

void usage(const std::string& bin_name) {
    std::cerr << "Usage : " << bin_name << " <connection_url>\n"
              << "Connection URL format should be :\n"
              << " For TCP : tcp://[server_host][:server_port]\n"
              << " For UDP : udp://[bind_host][:bind_port]\n"
              << " For Serial : serial:///path/to/serial/dev[:baudrate]\n"
              << "For example, to connect to the simulator use URL: udp://:14540\n";
}

std::shared_ptr<System> get_system(Mavsdk& mavsdk) {
    // std::cout << "Waiting to discover system...\n";
    auto prom = std::promise<std::shared_ptr<System>>{};
    auto fut = prom.get_future();

    // We wait for new systems to be discovered, once we find one that has an
    // autopilot, we decide to use it.
    mavsdk.subscribe_on_new_system([&mavsdk, &prom]() {
        auto system = mavsdk.systems().back();

        if (system->has_autopilot()) {
            // std::cout << "Discovered autopilot\n";

            // Unsubscribe again as we only want to find one system.
            mavsdk.subscribe_on_new_system(nullptr);
            prom.set_value(system);
        }
    });

    // We usually receive heartbeats at 1Hz, therefore we should find a
    // system after around 3 seconds max, surely.
    if (fut.wait_for(seconds(3)) == std::future_status::timeout) {
        std::cerr << "No autopilot found.\n";
        return {};
    }

    // Get discovered system now.
    return fut.get();
}

std::vector<float> create_6_6_mtx(std::vector<float>& tmp) {
    std::vector<float> pose_cov;

    for (int i = 0; i < 6; i++) {
        pose_cov.push_back(tmp[i]);
    }

    pose_cov.push_back(tmp[1]);

    for (int i = 6; i < 11; i++) {
        pose_cov.push_back(tmp[i]);
    }

    pose_cov.push_back(tmp[2]);
    pose_cov.push_back(tmp[7]);

    for (int i = 11; i < 15; i++) {
        pose_cov.push_back(tmp[i]);
    }

    pose_cov.push_back(tmp[3]);
    pose_cov.push_back(tmp[8]);
    pose_cov.push_back(tmp[12]);

    for (int i = 15; i < 18; i++) {
        pose_cov.push_back(tmp[i]);
    }

    pose_cov.push_back(tmp[4]);
    pose_cov.push_back(tmp[9]);
    pose_cov.push_back(tmp[13]);
    pose_cov.push_back(tmp[16]);

    for (int i = 18; i < 20; i++) {
        pose_cov.push_back(tmp[i]);
    }

    pose_cov.push_back(tmp[5]);
    pose_cov.push_back(tmp[10]);
    pose_cov.push_back(tmp[14]);
    pose_cov.push_back(tmp[17]);
    pose_cov.push_back(tmp[19]);

    pose_cov.push_back(tmp[20]);

    return pose_cov;
}

std::vector<float> create_3_3_mtx(std::vector<float>& tmp) {
    std::vector<float> ori_cov;

    ori_cov.push_back(tmp[15]); ori_cov.push_back(tmp[16]); ori_cov.push_back(tmp[17]);
    ori_cov.push_back(tmp[16]); ori_cov.push_back(tmp[18]); ori_cov.push_back(tmp[19]);
    ori_cov.push_back(tmp[17]); ori_cov.push_back(tmp[19]); ori_cov.push_back(tmp[20]);

    return ori_cov;
}


void calibrate_accelerometer(Calibration& calibration) {
    std::cout << "Calibrating accelerometer...\n";

    std::promise<void> calibration_promise;
    auto calibration_future = calibration_promise.get_future();

    calibration.calibrate_accelerometer_async(create_calibration_callback(calibration_promise));

    calibration_future.wait();
}

std::function<void(Calibration::Result, Calibration::ProgressData)>
create_calibration_callback(std::promise<void>& calibration_promise) {
    return [&calibration_promise](
               const Calibration::Result result, const Calibration::ProgressData progress_data) {
        switch (result) {
            case Calibration::Result::Success:
                std::cout << "--- Calibration succeeded!\n";
                calibration_promise.set_value();
                break;
            case Calibration::Result::Next:
                if (progress_data.has_progress) {
                    std::cout << "    Progress: " << progress_data.progress << '\n';
                }
                if (progress_data.has_status_text) {
                    std::cout << "    Instruction: " << progress_data.status_text << '\n';
                }
                break;
            default:
                std::cout << "--- Calibration failed with message: " << result << '\n';
                calibration_promise.set_value();
                break;
        }
    };
}

void calibrate_gyro(Calibration& calibration) {
    std::cout << "Calibrating gyro...\n";

    std::promise<void> calibration_promise;
    auto calibration_future = calibration_promise.get_future();

    calibration.calibrate_gyro_async(create_calibration_callback(calibration_promise));

    calibration_future.wait();
}

void calibrate_magnetometer(Calibration& calibration) {
    std::cout << "Calibrating magnetometer...\n";

    std::promise<void> calibration_promise;
    auto calibration_future = calibration_promise.get_future();

    calibration.calibrate_magnetometer_async(create_calibration_callback(calibration_promise));

    calibration_future.wait();
}

int main(int argc, char** argv) {
    if (argc != 2) {
        usage(argv[0]);
        return 1;
    }

    Mavsdk mavsdk;
    ConnectionResult connection_result = mavsdk.add_any_connection(argv[1]);

    if (connection_result != ConnectionResult::Success) {
        std::cerr << "Connection failed: " << connection_result << '\n';
        return 1;
    }

    auto system = get_system(mavsdk);
    if (!system) {
        return 1;
    }


    // auto calibration = Calibration(system);

    // calibrate_accelerometer(calibration);
    // calibrate_gyro(calibration);
    // calibrate_magnetometer(calibration);

    auto action = Action{system};
    auto mission = Mission{system};
    auto telemetry = Telemetry{system};

    while (!telemetry.health_all_ok()) {
        sleep_for(seconds(1));
    }


    const Action::Result arm_result = action.arm();

    if (arm_result != Action::Result::Success) {
        std::cerr << "Arming failed: " << arm_result << '\n';
        return 1;
    }

    Telemetry::Position home = telemetry.position();


    std::vector<Mission::MissionItem> mission_items;

    mission_items.push_back(make_mission_item(
        47.398170327054473,
        8.5456490218639658,
        10.0f,
        5.0f,
        false,
        20.0f,
        60.0f,
        Mission::MissionItem::CameraAction::None));

    mission_items.push_back(make_mission_item(
        47.398241338125118,
        8.5455360114574432,
        10.0f,
        2.0f,
        true,
        0.0f,
        -60.0f,
        Mission::MissionItem::CameraAction::TakePhoto));

    mission_items.push_back(make_mission_item(
        47.398139363821485,
        8.5453846156597137,
        10.0f,
        5.0f,
        true,
        -45.0f,
        0.0f,
        Mission::MissionItem::CameraAction::StartVideo));

    mission_items.push_back(make_mission_item(
        47.398058617228855,
        8.5454618036746979,
        10.0f,
        2.0f,
        false,
        -90.0f,
        30.0f,
        Mission::MissionItem::CameraAction::StopVideo));

    mission_items.push_back(make_mission_item(
        47.398100366082858,
        8.5456969141960144,
        10.0f,
        5.0f,
        false,
        -45.0f,
        -30.0f,
        Mission::MissionItem::CameraAction::StartPhotoInterval));

    // mission_items.push_back(make_mission_item(
    //     47.398001890458097,
    //     8.5455576181411743,
    //     10.0f,
    //     5.0f,
    //     false,
    //     0.0f,
    //     0.0f,
    //     Mission::MissionItem::CameraAction::StopPhotoInterval));

    mission_items.push_back(make_mission_item(
        home.latitude_deg,
        home.longitude_deg,
        home.relative_altitude_m,
        5.0f,
        false,
        0.0f,
        0.0f,
        Mission::MissionItem::CameraAction::TakePhoto));


    Mission::MissionPlan mission_plan{};
    mission_plan.mission_items = mission_items;
    const Mission::Result upload_result = mission.upload_mission(mission_plan);

    if (upload_result != Mission::Result::Success) {
        std::cerr << "Mission upload failed: " << upload_result << ", exiting.\n";
        return 1;
    }



    Mission::Result start_mission_result = mission.start_mission();


    /*==========================================================================

                                    GET TELEMETRY

    ==========================================================================*/
    std::ofstream fout("telemetry_new.txt");

    Telemetry::Odometry data;

    while (!mission.is_mission_finished().second) {
        data = telemetry.odometry();

        auto pose_cov = create_6_6_mtx(data.pose_covariance.covariance_matrix);
        auto ori_cov = create_3_3_mtx(data.pose_covariance.covariance_matrix);
        auto twist_cov = create_6_6_mtx(data.velocity_covariance.covariance_matrix);

        if (!isnan(pose_cov[0]) && !isnan(ori_cov[0]) && !isnan(twist_cov[0]) && !isnan(twist_cov[21])) {
            fout << data.q.w << " " << data.q.x << " " << data.q.y << " " << data.q.z << " ";
            fout << data.velocity_body.x_m_s << " " << data.velocity_body.y_m_s << " " << data.velocity_body.z_m_s << " ";
            fout << data.angular_velocity_body.roll_rad_s << " " << data.angular_velocity_body.pitch_rad_s << " " << data.angular_velocity_body.yaw_rad_s << " ";

            for (auto x : pose_cov) fout << x << " ";

            for (auto x : ori_cov) fout << x << " ";

            for (auto x : twist_cov) fout << x << " ";

            fout << data.time_usec << std::endl;
        }
    }


    const Action::Result hold_result = action.hold();
    if (hold_result != Action::Result::Success) {
        std::cout << "Failed to command Hold: " << hold_result << '\n';
        return 1;
    }

    const Action::Result land_result = action.land();
    if (land_result != Action::Result::Success) {
        std::cout << "Failed to command Land: " << land_result << '\n';
        return 1;
    }

    while (telemetry.armed()) {
        sleep_for(seconds(1));
    }

    const Action::Result disarm_result = action.disarm();
    if (disarm_result != Action::Result::Success) {
        std::cout << "Failed to command Disarm: " << disarm_result << '\n';
        return 1;
    }

    std::cout << "Disarmed" << std::endl;

    fout.close();
}
