import numpy as np
import matplotlib.pyplot as plt

# Returns an array with 1 in the region that contains the N centre-most grid points and 0 on the others.
def find_quasi_circle(N, grid_width, grid_height):

    current_grid = np.zeros((grid_height, grid_width))
    less_than_ideal_radius = np.sqrt(N/np.pi) / 1.2   
    n_cells_current = 0

    # Center of the grid in 0 based indexing
    a = grid_width/2 - 0.5
    b = grid_height/2 - 0.5

    # Increases the radius of the circle until if finds the correct size
    for radius in np.arange(less_than_ideal_radius, (min(grid_width, grid_height)/2), 0.1):
        last_grid = current_grid.copy()
        n_cells_last = n_cells_current

        # fill the circle
        for x in range(grid_width):
            for y in range(grid_height):
                if (x-a)**2 + (y-b)**2 < radius**2:
                    current_grid[x,y] = 1

        n_cells_current = np.count_nonzero(current_grid == 1)

        if n_cells_current == N:
            return current_grid
        
        if n_cells_current > N:

            n_cells_to_add = N - n_cells_last
            circle_border = np.ma.masked_where(last_grid == 1, current_grid)

            # Fill the outer cells in a diametrically opposit counter-clocksise way
            coordinates = np.where(circle_border == 1)
            sign = 1
            for i in range(n_cells_to_add):
                coord_to_add = [coordinates[0][i*sign], coordinates[1][i*sign]]
                last_grid[coord_to_add[0]][coord_to_add[1]] = 1
                sign *= -1
            
            return last_grid

# A test to see a 2D colormap of a given number of cells in the 201x201 grid
# width = 201
# height = 201
# resulting_quasi_circle = find_quasi_circle(97, width, height)
# plt.imshow(resulting_quasi_circle, interpolation='none')
# plt.show()

        
