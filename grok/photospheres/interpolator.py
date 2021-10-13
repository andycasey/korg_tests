
import numpy as np
import pickle
from scipy import interpolate
from tqdm import tqdm
import warnings

from .photosphere import Photosphere

class PhotosphereInterpolator(object):

    def __init__(self, photospheres, grid_keywords=None, interpolate_log_quantities=None):

        self.photospheres = photospheres

        # Build the grid of points.
        if grid_keywords is None:
            grid_keywords = tuple(self.photospheres[0].meta["grid_keywords"])
        self.grid_keywords = grid_keywords
        self.interpolate_log_quantities = interpolate_log_quantities or ()

        return None
        

    @property
    def grid_points(self):
        try:
            return self._grid_points
        except AttributeError:
            self._grid_points = np.array([
                [p.meta[k] for k in self.grid_keywords] for p in self.photospheres
            ])

            # Check for duplicates.
            remove_indices = []
            for i, column in enumerate(self._grid_points.T):
                if np.unique(column).size == 1:
                    # Warn, then remove.
                    warnings.warn(
                        f"Column index {i} ({self.grid_keywords[i]}) only has a single value: {column[0]}. "
                        f"Excluding it from interpolator dimensions."
                    )
                    remove_indices.append(i)
            
            if remove_indices:
                self.grid_keywords = tuple([kw for i, kw in enumerate(self.grid_keywords) if i not in remove_indices])
                N, D = self._grid_points.shape
                mask = np.ones(D, dtype=bool)
                mask[remove_indices] = False
                self._grid_points = self._grid_points[:, mask]
            
            unique = np.unique(self._grid_points, axis=0)
            if self._grid_points.shape != unique.shape:

                # Get an example.
                p_ = self._grid_points.view([('', self._grid_points.dtype)] * self._grid_points.shape[1])
                u_ = unique.view([('', unique.dtype)] * unique.shape[1])

                for each in u_:
                    match = (each == p_)
                    if sum(match) > 1:
                        example = f"Indices {tuple(np.where(match)[0])} have the same grid parameters: {each[0]}."
                        break
                        
                raise ValueError(
                    "There are duplicate points specified. It's likely that the photospheres have "
                    "additional values in the library that are not accounted for in the `grid_keywords` "
                    "meta. For example: the library of photospheres includes (teff, logg, fe_h, alpha_fe) "
                    "for each photosphere, but in the meta `grid_keywords` for each photosphere it is only"
                    "returning (teff, logg, fe_h) so there are multiple photospheres at each point. "
                    "For example:\n\n" + example + "\n\n "
                    "You can override the default `grid_keywords` when initiating the `PhotosphereInterpolator.`"
                )

        return self._grid_points


    def write(self, path):
        """
        Write this photosphere interpolator to disk. This will store all models and their metadata.

        :param path:
            The path to store the interpolator.
        """

        # Photosphere information first.
        column_names = self.photospheres[0].dtype.names

        N = len(self.photospheres)
        C = len(column_names)
        D = len(self.photospheres[0][column_names[0]])

        structure = np.nan * np.ones((N, C, D))
        for i in range(N):
            for j, column_name in enumerate(column_names):
                # In rare circumstances, a model can have one fewer depth points than it's neighbours.
                _structure = self.photospheres[i][column_name]
                structure[i, j, :len(_structure)] = _structure
                
        contents = {
            "points": self.grid_points,
            "column_names": column_names,
            "structure": structure,
            "meta": [p.meta for p in self.photospheres],
        }
        raise NotImplementedError()

        with open(path, "wb") as fp:
            pickle.dump(contents, fp)
        
        return None


    @classmethod
    def read(cls, path):
        """
        Read a photosphere interpolator from disk.
        
        :param path:
            The path where the library of models is stored.
        """

        with open(path, "rb") as fp:
            contents = pickle.load(fp)

        interpolator = cls([], None)
        interpolator._points = contents["points"]
        
        raise a



    def nearest_neighbour_indices(self, point, n):
        """
        Return the indices of the n nearest neighbours to the point.
        """

        distances = np.sum(((point - self.grid_points) / np.ptp(self.grid_points, axis=0))**2, axis=1)
        return distances.argsort()[:n]


    def __call__(self, method="linear", rescale=True, neighbours=30, **point):
        """
        Interpolate a photospheric structure at the given stellar parameters.
        """

        grid_points = self.grid_points

        missing_keys = set(self.grid_keywords).difference(point)
        if missing_keys:
            raise ValueError(f"Missing keyword arguments: {', '.join(missing_keys)}")
        
        opacity_column_name = "RHOX"
        interpolate_column_names = [n for n in self.photospheres[0].dtype.names[1:] if n != opacity_column_name]

        xi = np.array([point[k] for k in self.grid_keywords])

        lower, upper = (np.min(self.grid_points, axis=0), np.max(self.grid_points, axis=0))
        is_lower, is_upper = (xi < lower, xi > upper)
        if np.any(is_lower) or np.any(is_upper):
            is_bad = is_lower + is_upper
            indices = np.where(is_bad)[0]
            bad_values = xi[indices]
            raise ValueError(
                f"Point is outside the boundaries: {bad_values} (indices {indices}) outside bounds "
                f"(lower: {lower[indices]}, upper: {upper[indices]})"
            )

        grid_index = np.all(self.grid_points == xi, axis=1)
        if np.any(grid_index):
            grid_index = np.where(grid_index)[0][0]
            return self.photospheres[grid_index]

        # Work out what the optical depth points will be in our (to-be)-interpolated photosphere.
        neighbour_indices = self.nearest_neighbour_indices(xi, neighbours)

        # Protect Qhull from columns with a single value.
        cols = _protect_qhull(self.grid_points[neighbour_indices])  
        
        kwds = {
            "xi": xi[cols].reshape(1, len(cols)),
            "points": self.grid_points[neighbour_indices][:, cols],
            "values": np.array([self.photospheres[ni][opacity_column_name] for ni in neighbour_indices]),
            "method": method,
            "rescale": rescale
        }
        common_opacity_scale = interpolate.griddata(**kwds)

        # At the neighbouring N points, create splines of all the values
        # with respect to their own opacity scales, then calcualte the 
        # photospheric quantities on the common opacity scale.

        C = len(interpolate_column_names)
        N, D = kwds["values"].shape

        # Re-sample our neighbours onto the common opacity scale.
        neighbour_quantities = np.zeros((N, C, D))
        for i, neighbour_index in enumerate(neighbour_indices):
            neighbour = self.photospheres[neighbour_index]
            for j, column_name in enumerate(interpolate_column_names):
                tk = interpolate.splrep(neighbour[opacity_column_name], neighbour[column_name])
                neighbour_quantities[i, j, :] = interpolate.splev(common_opacity_scale, tk)[0]

        # TODO: interpoalte in log for some params
        kwds.update(values=neighbour_quantities)
        interpolated_quantities = interpolate.griddata(**kwds).reshape((C, D))

        # Get meta from neighbours.
        meta = self.photospheres[neighbour_indices[0]].meta.copy()

        # create a photosphere from these data.
        photosphere = Photosphere(
            data=np.vstack([common_opacity_scale, interpolated_quantities]).T,
            names=tuple([opacity_column_name] + interpolate_column_names),
            meta=meta
        )

        # Update the meta keywords for interpolated properties.
        for key in self.grid_keywords:
            photosphere.meta[key] = point[key]

        return photosphere


def _protect_qhull(a):
    return np.where([np.unique(a[:, i]).size > 1 for i in range(a.shape[1])])[0]
    