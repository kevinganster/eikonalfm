from tests.mesh import Mesh2D
import numpy as np
import eikonalfm
import matplotlib.pyplot as plt
import time

print("numpy version:", np.version.version)


def test_single(area, n, c_func, tau_exact, x_s=(0, 0), order=2):
    mesh = Mesh2D(*area, nx1=n)
    print(mesh)
    
    # source position
    x_s = mesh.pos_to_index(np.array(x_s))
    # velocity
    c = c_func(mesh.X)
    #np.savetxt("test.txt", c, fmt="%.2f", delimiter=", ", newline=",\n")
    
    tau = tau_exact(mesh.X, mesh.index_to_pos(x_s))
    print("exact solution:")
    print(np.array2string(tau.T, separator="    "))
    tau_fm = eikonalfm.fast_marching(c, x_s, mesh.dx, order)
    tau_ffm = eikonalfm.factored_fast_marching(c, x_s, mesh.dx, order)
    print("finished solving")
    print("ffm solution:")
    print(np.array2string(tau_ffm.T, separator="    "))

    print("error fast marching")
    print("inf-norm", np.linalg.norm(tau_fm - tau, ord=np.inf))
    print("2-norm", np.linalg.norm(tau_fm - tau, ord=2))
    print()
    print("error factored fast marching")
    print("inf-norm", np.linalg.norm(tau_ffm - tau, ord=np.inf))
    print("2-norm", np.linalg.norm(tau_ffm - tau, ord=2))
    
    
    levels = 10
    plt.figure()
    plt.suptitle("Fast Marching")
    ax = plt.subplot(3, 2, 1)
    ax.invert_yaxis()
    ax.set_title("Analytical solution")
    plot = ax.contourf(*mesh.X, tau, levels=levels)
#     plot = ax.pcolormesh(*mesh.X, tau)
    plt.colorbar(plot, ax=ax)
     
    ax = plt.subplot(3, 2, 2)
    ax.invert_yaxis()
    ax.set_title("Velocity c")
    plot = ax.contourf(*mesh.X, c, levels=levels)
    plt.colorbar(plot, ax=ax)
    
      
    ax = plt.subplot(3, 2, 3)
    ax.invert_yaxis()
    ax.set_title("Fast Marching")
    plot = ax.contourf(*mesh.X, tau_fm, levels=levels)
#     plot = ax.pcolormesh(*mesh.X, tau_fm)
    plt.colorbar(plot, ax=ax)
    # error
    ax = plt.subplot(3, 2, 4)
    ax.invert_yaxis()
    ax.set_title("Error")
    plot = ax.pcolormesh(*mesh.X, tau-tau_fm, cmap="seismic")
    plt.colorbar(plot, ax=ax)
       
    ax = plt.subplot(3, 2, 5)
    ax.invert_yaxis()
    ax.set_title("Factored Fast Marching")
    plot = ax.contourf(*mesh.X, tau_ffm, levels=levels)
#     plot = ax.pcolormesh(*mesh.X, tau_ffm)
    plt.colorbar(plot, ax=ax)
       
    ax = plt.subplot(3, 2, 6)
    ax.invert_yaxis()
    ax.set_title("Error")
    plot = ax.pcolormesh(*mesh.X, tau-tau_ffm, cmap="seismic")
    plt.colorbar(plot, ax=ax)
      
    plt.show()

    # phi = np.ones_like(c)
    # phi[tuple(x_s)] = 0
    # tau_skfmm = skfmm.travel_time(phi, c, dx=mesh.dx, order=2)


def test_isochrones(area, n, c_func, tau_exact, x_s=(0, 0), x_r=(1, 0), order=2):
    mesh = Mesh2D(*area, nx1=n)
    print(mesh)

    # velocity
    c = c_func(mesh.X)
    
    # source and receiver position
    x_s = mesh.pos_to_index(np.array(x_s))
    x_rec = mesh.pos_to_index(np.array(x_r))


    tau = tau_exact(mesh.X, mesh.index_to_pos(x_s))
    tau_fm = eikonalfm.fast_marching(c, x_s, mesh.dx, order)
    tau_ffm = eikonalfm.factored_fast_marching(c, x_s, mesh.dx, order)

    tau_rec = tau_exact(mesh.X, mesh.index_to_pos(x_rec))
    tau_rec_fm = eikonalfm.fast_marching(c, x_rec, mesh.dx, order)
    tau_rec_ffm = eikonalfm.factored_fast_marching(c, x_rec, mesh.dx, order)
    
    
    plt.figure()
    plt.suptitle("Isochrones")
    
    ax = plt.subplot(3, 2, 1)
    ax.invert_yaxis()
    ax.set_title("Exact solution")
    plot = ax.contour(*mesh.X, tau + tau_rec)
    ax.plot(*mesh.index_to_pos(x_s), color="red", marker="o", markersize=6)
    ax.plot(*mesh.index_to_pos(x_rec), color="red", marker="o", markersize=6)
    ax.clabel(plot, plot.levels, inline_spacing=0.1, fontsize=10)
    
    ax = plt.subplot(3, 2, 3)
    ax.invert_yaxis()
    ax.set_title("Fast Marching")
    plot = ax.contour(*mesh.X, tau_fm + tau_rec_fm)
    ax.plot(*mesh.index_to_pos(x_s), color="red", marker="o", markersize=6)
    ax.plot(*mesh.index_to_pos(x_rec), color="red", marker="o", markersize=6)
    ax.clabel(plot, plot.levels, fontsize=10)
    
    ax = plt.subplot(3, 2, 5)
    ax.invert_yaxis()
    ax.set_title("Factored Fast Marching")
    plot = ax.contour(*mesh.X, tau_ffm + tau_rec_ffm)
    ax.plot(*mesh.index_to_pos(x_s), color="red", marker="o", markersize=6)
    ax.plot(*mesh.index_to_pos(x_rec), color="red", marker="o", markersize=6)
    ax.clabel(plot, plot.levels, fontsize=10)
    
    plt.show()


def test_convergence(area, n_list, c_func, tau_exact, x_s=(0, 0), order=2):
    data = {
        'fm': {
            'inf': [],
            '2': [],
            'time' : []
        },
        'ffm': {
            'inf': [],
            '2': [],
            'time': []
        }
    }
    h_list = []

    for n in n_list:
        mesh = Mesh2D(*area, nx1=n)
        print("-" * 50)
        print(mesh)
        x = mesh.pos_to_index(np.array(x_s))
        h_list.append(mesh.dx[0])

        # velocity
        c = c_func(mesh.X)

        tau = tau_exact(mesh.X, mesh.index_to_pos(x))

        start = time.perf_counter()
        tau_fm = eikonalfm.fast_marching(c, x, mesh.dx, order)
        err = (tau - tau_fm).flatten()
        data['fm']['time'].append(time.perf_counter() - start)
        data['fm']['inf'].append(np.linalg.norm(err, ord=np.inf))
        data['fm']['2'].append(mesh.dx[0] * np.linalg.norm(err, ord=2))

        start = time.perf_counter()
        tau_ffm = eikonalfm.factored_fast_marching(c, x, mesh.dx, order)
        err = (tau - tau_ffm).flatten()
        data['ffm']['time'].append(time.perf_counter() - start)
        data['ffm']['inf'].append(np.linalg.norm(err, ord=np.inf))
        data['ffm']['2'].append(mesh.dx[0] * np.linalg.norm(err, ord=2))

        print()
        print("""Statistic:             fm     \t  ffm
        l_inf error: {:.2e},\t{:.2e}
        l_2 error:   {:.2e},\t{:.2e}
        time:        {:.2f}s,\t{:.2f}s""".format(
            data['fm']['inf'][-1], data['ffm']['inf'][-1],
            data['fm']['2'][-1], data['ffm']['2'][-1],
            data['fm']['time'][-1], data['ffm']['time'][-1]
        ))

    plt.figure()
    ax = plt.subplot(2, 3, 1)
    ax.set_title("Fast Marching")
    ax.set_xlabel("h")
    ax.plot(h_list, data['fm']['inf'], label="inf-norm")
    ax.plot(h_list, data['fm']['2'], label="2-norm")
    ax.set_xticks(h_list)
    ax.legend()

    ax = plt.subplot(2, 3, 2)
    ax.set_title("Factored Fast Marching")
    ax.set_xlabel("h")
    ax.plot(h_list, data['ffm']['inf'], label="inf-norm")
    ax.plot(h_list, data['ffm']['2'], label="2-norm")
    ax.set_xticks(h_list)
    ax.legend()

    ax = plt.subplot(2, 3, 3)
    ax.set_title("Time")
    ax.set_xlabel("n")
    ax.set_ylabel("sec")
    ax.loglog(n_list, data['fm']['time'], label="fast marching")
    ax.loglog(n_list, data['ffm']['time'], label="factored fast marching")
#     ax.loglog(n_list, [1e-4 * n * np.log(n) for n in n_list], label="n log(n)")
    ax.legend()
    
    
    # logplots
    ax = plt.subplot(2, 3, 4)
    ax.set_title("Fast Marching")
    ax.set_xlabel("h")
    ax.loglog(h_list, data['fm']['inf'], label="inf-norm")
    ax.loglog(h_list, data['fm']['2'], label="2-norm")
    ax.legend()

    ax = plt.subplot(2, 3, 5)
    ax.set_title("Factored Fast Marching")
    ax.set_xlabel("h")
    ax.loglog(h_list, data['ffm']['inf'], label="inf-norm")
    ax.loglog(h_list, data['ffm']['2'], label="2-norm")
    ax.legend()
    
    plt.show()